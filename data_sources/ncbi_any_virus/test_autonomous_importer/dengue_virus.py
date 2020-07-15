import os
import re
import time
from collections import Counter, OrderedDict
from concurrent.futures.thread import ThreadPoolExecutor
from datetime import datetime
from decimal import Decimal
from typing import Callable, List, Optional, Tuple

import lxml
from Bio import Entrez
from dateutil.parser import parse
from tqdm import tqdm

import database_tom
import vcm
import xml_helper
from geo_groups import geo_groups
from xml_helper import text_at_node

Entrez.email = "Your.Name.Here@example.org"
from loguru import logger
from lxml import etree

from locations import get_local_folder_for, FileType


reference_sequences = {}    # a dict organism_name/alias -> (reference sample, organism taxon)
organism_aliases = {}       # map principal_organism_name -> aliases
virus_ids = {}              # map of organism_name/alias -> virus_id
DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = 30


class AnyNCBITaxon:

    def __init__(self, xml_tree_file_path: str):
        self.tax_tree: etree.ElementTree = \
            etree.parse(xml_tree_file_path, parser=etree.XMLParser(remove_blank_text=True)) \
            .xpath('/TaxaSet/Taxon')[0]

    def __str__(self):
        return f'{self.taxon_id()}.xml'

    def taxon_id(self):
        return text_at_node(self.tax_tree, './TaxId', mandatory=True)

    def taxon_name(self):
        return text_at_node(self.tax_tree, './ScientificName')

    def family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName')

    def sub_family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName', mandatory=False)

    def genus(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName')

    def species(self):
        # species_taxon_id = text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/TaxId')
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName', mandatory=False)

    def equivalent_names(self):
        genbank_acronym = text_at_node(self.tax_tree, './/GenbankAcronym', mandatory=False)
        equivalent_names = self.tax_tree.xpath('.//EquivalentName')
        equivalent_names = [x.text for x in equivalent_names]
        if genbank_acronym:
            equivalent_names.insert(0, genbank_acronym)
        equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
        equivalent_names = ", ".join(equivalent_names)
        return equivalent_names


class AnyNCBIVNucSample:
    """
    Set of getters to take the relevant information from a virus sample in INSDSeq XML format.
    Example of virus sample @ https://www.ncbi.nlm.nih.gov/nuccore/MN908947
    INSDSeq XML format DTD @ https://www.ncbi.nlm.nih.gov/data_specs/dtd/INSD_INSDSeq.mod.dtd
    """

    re_structured_comment = re.compile(r'##Assembly-Data-START## ;(.*); ##Assembly-Data-END##')
    cached_taxon_id = {}
    re_anna_RuL3z = re.compile(r'(\d+)([,|.]?)(\d*).*')
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)
    re_patent_submission_date_lab = re.compile(r'^.* (\d+-\w+-\d+)(?: (.+))?$')

    def __init__(self, virus_sample_file_path: str, internal_accession_id):
        """
        :param virus_sample_file_path: virus sample file (INDSeq XML) path
        :param internal_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
        Raises lxml.etree.XMLSyntaxError if the source file is malformed or empty.
        """
        super().__init__()
        self.sample_xml: lxml.etree.ElementTree = etree.parse(source=virus_sample_file_path,
                                                              parser=etree.XMLParser(remove_blank_text=True))
        self._internal_accession_id = internal_accession_id

        # init internal variables
        self._journal = None
        self._host_value_list = None
        self._annotations = None
        self._nuc_seq = None
        self.taxonomy_info: Optional[AnyNCBITaxon] = None

    def __str__(self):
        return f'{self.internal_id()}.xml'

    def set_taxonomy_info(self, taxonomy_info: AnyNCBITaxon):
        self.taxonomy_info = taxonomy_info

    def __getattr__(self, name):
        """This is invoked after a lookup of {name} into this class. This provides a callback for a fallback action other
        than raising AttributeError"""
        return getattr(self.taxonomy_info, name)

    def internal_id(self):
        """
        An id used internally by this program when needed to log some errors. Choose whatever form and type you like.
        This id is not going to be written in the output data and possibly it will be printed only on the console in order
        to trace errors.
        """
        return self._internal_accession_id

    def primary_accession_number(self):
        return text_at_node(self.sample_xml, './/INSDSeq_accession-version')

    def alternative_accession_number(self):
        return str(self._internal_accession_id) if self._internal_accession_id else None

    def organism(self):
        return text_at_node(self.sample_xml, './/INSDSeq_organism', True)

    def strain(self):
        strain = text_at_node(self.sample_xml,
                              './/INSDQualifier[./INSDQualifier_name/text() = "strain"]/INSDQualifier_value', False)
        if not strain:
            # try to get the strain from the filed FEATURES -> SOURCE -> /isolate
            strain = text_at_node(self.sample_xml,
                                  './/INSDQualifier[./INSDQualifier_name/text() = "isolate"]/INSDQualifier_value',
                                  False)
        return strain

    def is_reference(self):
        return text_at_node(self.sample_xml, './/INSDKeyword', mandatory=False) == 'RefSeq'

    def get_reference_sequence(self):
        try:
            return reference_sequences[self.organism()]
        except KeyError as e:
            logger.error(f'Error in sample {self.internal_id()}: organism {self.organism()} does not have a matching reference sequence.'
                         f'reference sequences are available for organisms {list(reference_sequences.keys())}')
            raise e

    def is_complete(self):
        definition = text_at_node(self.sample_xml, ".//INSDSeq_definition")
        definition_0 = definition.split(";")[0]
        definition_0_last = definition_0.split(",")[-1]
        definition_0_last = definition_0_last.strip()
        if definition_0_last in ['complete genome', ]:
            return True
        elif definition_0_last in ['partial cds', 'complete cds', 'partial genome']:
            return False
        else:
            length = self.length()
            reference_sample, _ = self.get_reference_sequence()
            reference_length = reference_sample.length()
            if length is not None and reference_length is not None and length < 0.95*reference_length:
                return False
            else:
                return None

    def nucleotide_sequence(self):
        if not self._nuc_seq:
            try:
                self._nuc_seq = text_at_node(self.sample_xml, './/INSDSeq_sequence')
            except AssertionError:
                raise ValueError(f'Insertion of sequence {self.internal_id()} skipped because of missing nucleotide seq.')
        return self._nuc_seq

    @staticmethod
    def strand():
        return 'positive'

    def is_single_stranded(self):
        strandedness = text_at_node(self.sample_xml, './INSDSeq/INSDSeq_strandedness', mandatory=False)
        if strandedness and strandedness.lower() == 'single':
            return True
        else:
            logger.warning(f'unable to determine strandedness from string {strandedness} in sample {self.internal_id()}')
            return None

    def molecule_type(self):
        return text_at_node(self.sample_xml, './INSDSeq/INSDSeq_moltype', mandatory=False)

    @staticmethod
    def is_positive_stranded():
        return True

    def length(self):
        length = int(text_at_node(self.sample_xml, './/INSDSeq_length'))
        assert length == len(self.nucleotide_sequence())
        return length

    def gc_percent(self):
        c = Counter(self.nucleotide_sequence().lower())
        count_known_nucleotides = (c['g'] + c['c'] + c['a'] + c['t'] + c['u'])
        if count_known_nucleotides != 0:
            gc_percentage = (c['g'] + c['c']) / count_known_nucleotides * 100
            gc_percentage = Decimal(gc_percentage)
            gc_percentage = round(gc_percentage, 2)
            return gc_percentage
        else:
            return 0

    def n_percent(self):
        length = self.length()
        if length != 0:
            c = Counter(self.nucleotide_sequence().lower())
            n_percentage = (c['n']) / length * 100
            n_percentage = Decimal(n_percentage)
            n_percentage = round(n_percentage, 2)
            return n_percentage
        else:
            return 0

    @staticmethod
    def lineage():
        return None

    @staticmethod
    def clade():
        return None

    def sequencing_technology(self):
        return self._structured_comment(self.sample_xml, 'Sequencing Technology')

    def assembly_method(self):
        return self._structured_comment(self.sample_xml, 'Assembly Method')

    def coverage(self) -> Optional[str]:
        _input = self._structured_comment(self.sample_xml, 'Coverage')
        if not _input:
            return None
        info = self.re_anna_RuL3z.match(_input)
        if not info:
            raise ValueError(f'coverage string {_input} doesn\'t match the regex.')
        else:
            o1 = info.group(1)  # integer part
            o3 = info.group(3)  # decimal (optional)
            # apply Anna's rules
            output = int(o1)
            if o3 and len(o3) < 3:
                decimals = round(int(o3) / 10)
                output += decimals
            return str(output)

    def collection_date(self):
        collection_date = text_at_node(self.sample_xml,
                                       '..//INSDQualifier[./INSDQualifier_name/text() = "collection_date"]/INSDQualifier_value',
                                       mandatory=False)
        if collection_date:
            collection_date = parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
        return collection_date

    def submission_date(self):
        try:
            return datetime.strptime(self._init_and_get_journal()[0], '%d-%b-%Y')
        except TypeError:
            return None

    def isolation_source(self):
        source = text_at_node(self.sample_xml,
                              '..//INSDQualifier[./INSDQualifier_name/text() = "isolation_source"]/INSDQualifier_value',
                              mandatory=False)
        if source:
            source = source.lower()
        return source

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        country = None
        region = None
        node = text_at_node(self.sample_xml,
                            '..//INSDQualifier[./INSDQualifier_name/text() = "country"]/INSDQualifier_value',
                            mandatory=False)
        if node:
            node = node.split(":")
            country = node[0]
            region = node[1] if len(node) > 1 else None
        geo_group = geo_groups.get(country.lower()) if country else None
        return country, region, geo_group

    def _init_and_get_journal(self):
        if not self._journal:
            try:
                references = self.sample_xml.xpath('.//INSDReference[./INSDReference_title/text() = "Direct Submission"]')
                assert len(references) > 0, 'there must be at least one direct submission'
                reference = references[0]
                journal_string = text_at_node(reference, "./INSDReference_journal")
                assert journal_string.startswith(
                    "Submitted "), 'Cannot find submitted in the Journal of direct submission reference'
                self._journal = re.split("[()]", journal_string, maxsplit=2)[1:]
                assert len(self._journal) == 2, f"Journal value problem '{journal_string}' {self._journal}"

                #   # journal details
                #     authors = reference.xpath('.//INSDAuthor')
                #     authors = [x.text for x in authors]
                #     authors = ", ".join(authors)
                #     title = text_at_node(reference, "./INSDReference_title")
                #     journal = text_at_node(reference, "./INSDReference_journal")
                #     publication_date = None
                #     pubmed_id = text_at_node(reference, "./INSDReference_pubmed" , mandatory=False)
                #     popset = None
            except AssertionError as e:
                # try to get patent information
                definition = text_at_node(self.sample_xml, './/INSDSeq_definition', mandatory=False)
                if definition and 'patent' in definition.lower():
                    references = self.sample_xml.xpath('.//INSDReference/INSDReference_journal')
                    assert len(references) > 0, 'no journal info for patent'
                    for ref in references:
                        journal_string: str = text_at_node(ref, xpath_string='.', mandatory=True)
                        _journal_parts = self.re_patent_submission_date_lab.match(journal_string)
                        if _journal_parts:
                            self._journal = _journal_parts.groups()
                            break
        return self._journal

    @staticmethod
    def originating_lab() -> Optional[str]:
        return None

    def sequencing_lab(self) -> Optional[str]:
        try:
            return self._init_and_get_journal()[1]
        except TypeError:
            return None

    def _init_and_get_host_values(self):
        if not self._host_value_list:
            host_node = text_at_node(self.sample_xml,
                                     './/INSDQualifier[./INSDQualifier_name/text() = "host"]/INSDQualifier_value',
                                     mandatory=False)
            self._host_value_list = [x.strip() for x in host_node.split(";")] if host_node else None
            if self._host_value_list:
                for i in range(1, len(self._host_value_list)):
                    self._host_value_list[i] = self._host_value_list[i].lower()
        return self._host_value_list

    # host taxon name
    def host_taxon_name(self) -> Optional[str]:
        host = self._init_and_get_host_values()
        return host[0] if host else None

    # host taxon id
    def host_taxon_id(self) -> Optional[int]:
        taxon_name = self.host_taxon_name()
        if not taxon_name:
            return None
        else:
            taxon_id = self.cached_taxon_id.get(taxon_name)
            if taxon_id == -1:  # -1 means the cached taxon_id for this taxon name was searched before
                return None
            elif taxon_id is None:
                try:
                    with Entrez.esearch(db="taxonomy", term=f'"{taxon_name}"[Scientific Name]', rettype=None,
                                        retmode="xml") as handle:
                        response = Entrez.read(handle)
                        if response['Count'] == '1':
                            taxon_id = int(response['IdList'][0])
                            self.cached_taxon_id[taxon_name] = taxon_id
                        else:
                            logger.warning(f'can\'t find the taxon id for taxon name {taxon_name}')
                            self.cached_taxon_id[
                                taxon_name] = -1  # save -1 in cache to distinguish from non cached taxon_ids
                            taxon_id = None
                except:
                    logger.exception(
                        f'Exception occurred while fetching the taxon id of {taxon_name} in sample {self.internal_id()}')
            return taxon_id

    def gender(self) -> Optional[str]:
        host = self._init_and_get_host_values()
        if host:
            gender = 'male' if ('male' in host) else 'female' if 'female' in host else None
            for h in host:
                words = h.split()
                if 'm' in words:
                    gender = 'male'
                    break
                elif 'f' in words:
                    gender = 'female'
                    break
            return gender
        else:
            return None

    def age(self) -> Optional[int]:
        host = self._init_and_get_host_values()
        if host:
            age = next(filter(lambda x: 'age' in x, host), None)
            if age:
                age = int(next(filter(lambda x: x.isdigit(), age.split()), None))
            return age
        else:
            return None

    def database_source(self) -> str:
        return 'RefSeq' if self.is_reference() else 'GenBank'

    def bioproject_id(self):
        return text_at_node(self.sample_xml, './/INSDXref[./INSDXref_dbname/text() = "BioProject"]/INSDXref_id',
                            mandatory=False)

    @classmethod
    def _structured_comment(cls, el, key):
        comment = text_at_node(el, './/INSDSeq_comment', False)
        if comment:
            sub = cls.re_structured_comment.findall(comment)
            assert len(sub) <= 1, f"multiple structured_comment {len(sub)}"
            if sub:
                sub = sub[0]
                subs = sub.split(";")
                subs = [x for x in subs if key in x]
                assert len(subs) <= 1, f"multiple structured_comment for key: {key} {len(subs)}"
                if subs:
                    return subs[0].split("::")[1].strip()
                else:
                    return None
            else:
                return None
        else:
            return None


# noinspection PyPep8Naming
def reference_samples_from_organism(ncbi_taxon_id: int, log_with_name: str, TaxonomyWrapperClass=AnyNCBITaxon, SampleWrapperClass=AnyNCBIVNucSample):
    logger.info(f'importing reference samples and taxonomy of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of refseq samples...')
    refseq_accession_ids = _get_refseq_samples_accession_ids(ncbi_taxon_id)

    logger.trace(f'download and processing of reference sequences...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    download_taxonomy_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.TaxonomyData)
    for sample_id in tqdm(refseq_accession_ids):
        sample_path = _download_or_get_virus_sample_as_xml(download_sample_dir_path, sample_id)
        ref_sample = SampleWrapperClass(sample_path, sample_id)

        principal_organism_taxon_name = ref_sample.organism()
        logger.trace(f'downloading taxonomy data of organism {principal_organism_taxon_name}...')
        organism_taxon_id, related_organism_taxa_ids = _find_taxon_ids(principal_organism_taxon_name)
        taxonomy_file = _download_virus_taxonomy_as_xml(download_taxonomy_dir_path, organism_taxon_id)
        principal_organism = TaxonomyWrapperClass(taxonomy_file)

        # building a cache of reference sequences and organisms for all the samples
        reference_sequences[principal_organism_taxon_name] = (ref_sample, principal_organism)
        organism_aliases[principal_organism_taxon_name] = [principal_organism_taxon_name]
        for txid in related_organism_taxa_ids:
            related_org_taxonomy_file = _download_virus_taxonomy_as_xml(download_taxonomy_dir_path, txid)
            related_org = TaxonomyWrapperClass(related_org_taxonomy_file)
            reference_sequences[related_org.taxon_name()] = (ref_sample, principal_organism)
            organism_aliases[principal_organism_taxon_name].append(related_org.taxon_name())

        ref_sample.set_taxonomy_info(principal_organism)
        yield ref_sample


# noinspection PyPep8Naming
def other_samples_from_organism(ncbi_taxon_id: int, log_with_name: str, SampleWrapperClass=AnyNCBIVNucSample):
    logger.info(f'importing non-reference samples of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of non-refseq samples...')
    non_refseq_accession_ids = _get_other_samples_accession_ids(ncbi_taxon_id)
    # non_refseq_accession_ids = non_refseq_accession_ids[:10]  # TODO remove

    logger.trace(f'download and processing of non-reference sequences...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    for sample_id in tqdm(non_refseq_accession_ids):
        sample_path = _download_or_get_virus_sample_as_xml(download_sample_dir_path, sample_id)
        other_sample = SampleWrapperClass(sample_path, sample_id)
        yield other_sample


def _find_taxon_ids(organism_name: str):
    def do():
        # find the taxon id for this taxon organism
        with Entrez.esearch(db="taxonomy", term=f'"{organism_name}"[subtree]') as handle:
            response = Entrez.read(handle)
        principal_taxon_id = response['IdList'][-1]  # if multiple organism, take the last one (usually is the parent of all the organisms)
        sub_specie_taxa_ids = response['IdList'][:-1]   # related subtypes of the principal organism (can be an empty list)
        return principal_taxon_id, sub_specie_taxa_ids
    return _try_or_wait(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _download_virus_taxonomy_as_xml(containing_directory: str, taxon_id: int) -> str:
    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{os.path.sep}{taxon_id}.xml"
        if not os.path.exists(destination_file_path):
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                    open(destination_file_path, 'w') as f:
                for line in handle:
                    f.write(line)
        return destination_file_path
    return _try_or_wait(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _download_or_get_virus_sample_as_xml(containing_directory: str, sample_accession_id: int) -> str:
    """
    :param containing_directory: directory where the file will be downloaded and cached
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    def do():
        local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
        with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml") as handle, open(local_file_path, 'w') as f:
            for line in handle:
                f.write(line)
        return local_file_path
    return _try_or_wait(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _get_refseq_samples_accession_ids(taxon_id: int):
    def do():
        # total number of sequences
        with Entrez.esearch(db="nuccore",
                            term=f"(txid{taxon_id}[Organism:exp]) AND srcdb_refseq[Properties]",
                            rettype='count') as handle:
            response = Entrez.read(handle)
        total_records = int(response['Count'])
        # get accession ids
        with tqdm(total=total_records) as progress_bar:
            with Entrez.esearch(db="nuccore", term=f"(txid{taxon_id}[Organism:exp]) AND srcdb_refseq[Properties]") as handle:
                response = Entrez.read(handle)
                # # EXAMPLE RESPONSE
                # {
                #     "Count": "1",                       <-- number of total records matching the query
                #     "RetMax": "1",
                #     "RetStart": "0",
                #     "IdList": ["1798174254"],            <-- accession id of refseq
                #     "TranslationSet": [],
                #     "TranslationStack": [
                #         {"Term": "txid2697049[Organism]", "Field": "Organism", "Count": "5511", "Explode": "Y"},
                #         {"Term": "srcdb_refseq[Properties]", "Field": "Properties", "Count": "66995641", "Explode": "N"},
                #         "AND"
                #     ],
                #     "QueryTranslation": "txid2697049[Organism] AND srcdb_refseq[Properties]"
                # }
            refseq_acc_ids = [int(_id) for _id in response['IdList']]
            progress_bar.update(len(refseq_acc_ids))
            assert total_records == len(refseq_acc_ids), "no reference sample found. Please check: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch"
            return refseq_acc_ids
    return _try_or_wait(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _get_other_samples_accession_ids(taxon_id: int) -> (List[int], List[int]):
    def do():
        # DO PAGINATION
        # total number of sequences
        with Entrez.esearch(db="nuccore",
                            term=f"(txid{taxon_id}[Organism:exp]) NOT srcdb_refseq[Properties]",
                            rettype='count') as handle:
            response = Entrez.read(handle)
            total_records = int(response['Count'])
        # get pages
        non_refseq_accessions_ids = list()
        RECORDS_PER_PAGE = 1000
        page_number = 0
        with tqdm(total=total_records) as progress_bar:
            while total_records > page_number * RECORDS_PER_PAGE:
                with Entrez.esearch(db="nuccore",
                                    term=f"(txid{taxon_id}[Organism:exp]) NOT srcdb_refseq[Properties]",
                                    retmax=RECORDS_PER_PAGE, retstart=page_number * RECORDS_PER_PAGE) as handle:
                    response = Entrez.read(handle)
                for x in response['IdList']:
                    non_refseq_accessions_ids.append(int(x))
                page_number += 1
                progress_bar.update(
                    RECORDS_PER_PAGE if page_number * RECORDS_PER_PAGE < total_records else total_records - (
                                (page_number - 1) * RECORDS_PER_PAGE))
        if len(non_refseq_accessions_ids) != total_records:
            raise IOError('Some of the non-refseq accession ids were not correctly downloaded')
        return non_refseq_accessions_ids
    return _try_or_wait(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _try_or_wait(n_times: int, or_wait_secs: int, function: Callable, *args, **kwargs):
    with ThreadPoolExecutor(max_workers=1) as executor:
        try:
            return executor.submit(function, *args, **kwargs).result()
        except IOError as e:
            n_times -= 1
            if n_times > 0:
                logger.error(f'Error while invoking {function.__name__} with args: {args}\nkwargs: {kwargs}. '
                             f'New attempt in {or_wait_secs} secs.'
                             f'Reason of error: {str(type(e))} {e.args}')
                executor.submit(time.sleep, or_wait_secs).result()
                return _try_or_wait(n_times, or_wait_secs, function, *args, **kwargs)
            else:
                raise e


#   #   MAIN PIPELINE
def import_into_vcm_all_except_annotations_nuc_vars(session: database_tom.Session, sample: AnyNCBIVNucSample):
    virus_id = virus_ids[sample.organism()]  # the virus_id is the one associated with the reference sequence of the principal organism
    experiment = vcm.create_or_get_experiment(session, sample)
    host_sample = vcm.create_or_get_host_sample(session, sample)
    sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
    sequence = vcm.create_or_get_sequence(session, sample, virus_id, experiment, host_sample, sequencing_project)


def import_references_into_vcm(session: database_tom.Session, sample: AnyNCBIVNucSample):
    virus_id = vcm.create_or_get_virus(session, sample).virus_id

    # associate the organism of the reference sequence and related organisms to this virus_id
    this_org_aliases: [str] = organism_aliases[sample.organism()]
    logger.debug(f'organism associated to {sample.organism()} are {this_org_aliases}')
    global virus_ids
    for organism in this_org_aliases:
        virus_ids[organism] = virus_id
    logger.debug(f'virus_ids: {virus_ids}')

    import_into_vcm_all_except_annotations_nuc_vars(session, sample)


for refseq in reference_samples_from_organism(12637, 'Dengue virus'):
    logger.debug(organism_aliases)
    database_tom.try_py_function(
        import_references_into_vcm, refseq)

for other_seq in other_samples_from_organism(12637, 'Dengue virus'):
    try:
        database_tom.try_py_function(
            import_into_vcm_all_except_annotations_nuc_vars, other_seq
        )
    except:
        logger.exception(f'exception occurred while working on virus sample {other_seq}')



