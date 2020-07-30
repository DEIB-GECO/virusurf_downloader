import os
from os.path import sep
import re
import time
import urllib
from collections import Counter, OrderedDict
from concurrent.futures.thread import ThreadPoolExecutor
from queuable_jobs import Job, Boss, Worker
from datetime import datetime
from decimal import Decimal
from typing import Callable, List, Optional, Tuple

import dateutil.parser as dateparser
import lxml
from tqdm import tqdm
import database_tom
import vcm_mock as vcm
from geo_groups import geo_groups
from xml_helper import text_at_node
from loguru import logger
from lxml import etree
from locations import get_local_folder_for, FileType
from Bio import Entrez
import pickle
Entrez.email = "Your.Name.Here@example.org"
from pipeline_nuc_variants__annotations__aa import sequence_aligner

nucleotide_reference_sequence: Optional[str] = None # initialized elsewhere
annotation_file_path: Optional[str] = None  # initialized elsewhere
virus_sequence_chromosome_name: Optional[str] = None  # initialized elsewhere
snpeff_db_name: Optional[str] = None  # initialized elsewhere
log_with_name: Optional[str] = None
DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = 30


class AnyNCBITaxon:

    def __init__(self, xml_tree_file_path: str):
        self.tax_tree: etree.ElementTree = \
            etree.parse(xml_tree_file_path, parser=etree.XMLParser(remove_blank_text=True)) \
            .xpath('/TaxaSet/Taxon')[0]
        rank = self.tax_tree.xpath('./Rank')     # xpath returns a list also for single nodes
        if rank:
            self.suggested_from_other_method[text_at_node(rank[0], '.').lower()] = self.taxon_name()

    suggested_from_other_method = {}

    def __str__(self):
        return f'{self.taxon_id()}.xml'

    def taxon_id(self):
        return text_at_node(self.tax_tree, './TaxId', mandatory=True)

    def taxon_name(self):
        return text_at_node(self.tax_tree, './ScientificName')

    def family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName') \
               or self.suggested_from_other_method.get('family')

    def sub_family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName', mandatory=False) \
               or self.suggested_from_other_method.get('subfamily')

    def genus(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName') \
               or self.suggested_from_other_method.get('genus')

    def species(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName', mandatory=False) \
               or self.suggested_from_other_method.get('species')

    def equivalent_names(self):
        genbank_acronym = text_at_node(self.tax_tree, './/GenbankAcronym', mandatory=False)
        equivalent_names = self.tax_tree.xpath('.//EquivalentName')
        equivalent_names = [x.text for x in equivalent_names]
        if genbank_acronym:
            equivalent_names.insert(0, genbank_acronym)
        equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
        equivalent_names = ", ".join(equivalent_names)
        return equivalent_names

    @staticmethod
    def molecule_type():
        return 'RNA'

    @staticmethod
    def is_single_stranded():
        return True

    def is_positive_stranded(self):
        if 'ebolavirus' in self.taxon_name().lower():
            return False
        else:
            return True


class AnyNCBIVNucSample:
    """
    Set of getters to take the relevant information from a virus sample in INSDSeq XML format.
    Example of virus sample @ https://www.ncbi.nlm.nih.gov/nuccore/MN908947
    INSDSeq XML format DTD @ https://www.ncbi.nlm.nih.gov/data_specs/dtd/INSD_INSDSeq.mod.dtd
    """

    re_structured_comment = re.compile(r'##Assembly-Data-START## ;(.*); ##Assembly-Data-END##')
    cached_taxon_id = {}
    re_anna_RuL3z = re.compile(r'\D*(\d+)([,|.]?)(\d*).*')
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)
    re_patent_submission_date_lab = re.compile(r'^.* (\d+-\w+-\d+)(?: (.+))?$')
    gender_suggested_by_other_method: Optional[str] = None
    country_suggested_by_other_method: Optional[str] = None
    host_name_suggested_by_other_method:Optional[str] = None

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
            try:
                strain = text_at_node(self.sample_xml,
                                      './/INSDQualifier[./INSDQualifier_name/text() = "isolate"]/INSDQualifier_value',
                                      False)
            except AssertionError as e:
                logger.warning(f'In sample {self.internal_id()}: {str(e)}.')
                isolates = self.sample_xml.xpath('.//INSDQualifier[./INSDQualifier_name/text() = "isolate"]/INSDQualifier_value')
                fallback = text_at_node(isolates[0], '.', mandatory=True)
                logger.warning(f'Isolate strings are {isolates}. Only {fallback} was retained as strain')
                return fallback
        return strain

    def is_reference(self):
        keyword_nodes = self.sample_xml.xpath('.//INSDSeq_keywords/INSDKeyword')
        return 'RefSeq' in [text_at_node(node, '.', mandatory=False) for node in keyword_nodes]

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
            reference_length = len(nucleotide_reference_sequence)
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
            logger.warning(f'coverage string in sample {self.internal_id()} is {_input} and cannot be parsed. '
                           f'EXPERIMENT_TYPE.coverage set to NULL.\n'
                           f'Full comment content is: {text_at_node(self.sample_xml, ".//INSDSeq_comment", False)}')
            return None
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
            try:
                collection_date = dateparser.parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
            except dateparser._parser.ParserError as e:
                if '/' in collection_date:
                    fallback = dateparser.parse(collection_date.split('/')[0], default=self.default_datetime).strftime('%Y-%m-%d')
                else:
                    fallback = None
                logger.warning(f'collection date string"{collection_date}" was parsed as {fallback}')
                return fallback
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
            # check if it contains gender info
            if 'female' in source:
                self.gender_suggested_by_other_method = 'female'
                if source == 'female':
                    return None
            elif 'male' in source:
                self.gender_suggested_by_other_method = 'male'
                if source == 'male':
                    return None
            elif 'sapiens' in source:
                self.host_name_suggested_by_other_method = 'Homo sapiens'
                return None
            # check if it contains country info
            country = geo_groups.get(source.lower())
            if country:
                self.country_suggested_by_other_method = country
                return None

        return source

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        country = None
        region = None
        node = text_at_node(self.sample_xml,
                            '..//INSDQualifier[./INSDQualifier_name/text() = "country"]/INSDQualifier_value',
                            mandatory=False)
        if node:
            node = node.split(":")
            country = node[0] or self.country_suggested_by_other_method     # can be None
            region = node[1].strip() if len(node) > 1 else None
        else:
            country = self.country_suggested_by_other_method    # can be None
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
            except AssertionError:
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
        if not host:
            return None
        host = host[0]
        # remove stuff in parentheses
        try:
            first_par = host.index('(')     # raises ValueError if ( is not found
            after_second_par = host.index(')', first_par) + 1
            output = host[:first_par].strip()
            if after_second_par < len(host):
                output = output + ' ' + host[after_second_par:].strip()
            return output
        except ValueError:
            return host

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
                    with Entrez.esearch(db="taxonomy", term=taxon_name, rettype=None,
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
            return gender or self.gender_suggested_by_other_method      # can be None
        else:
            return self.gender_suggested_by_other_method    # can be None

    def age(self) -> Optional[int]:
        host = self._init_and_get_host_values()
        if host:
            age = next(filter(lambda x: 'age' in x, host), None)
            if age:
                age = int(''.join(filter(lambda x: x.isdigit(), age)))
                # age = int(next(filter(lambda x: x.isdigit(), age.split()), None))
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

    def on_before_multiprocessing(self):
        self.nucleotide_sequence()  # make sure nucleotide variants are cached in this sample
        self.sample_xml = None  # release xml reference to prevent errors with multiprocessing


# noinspection PyPep8Naming
def _all_samples_from_organism(samples_query: str, log_with_name: str, SampleWrapperClass=AnyNCBIVNucSample, from_sample: Optional[int] = None, to_sample: Optional[int] = None):
    logger.trace(f'importing samples of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of samples...')
    accession_ids = _get_samples_accession_ids(samples_query)
    if from_sample is not None and to_sample is not None:
        accession_ids.sort()
        accession_ids = accession_ids[from_sample:to_sample]

    logger.trace(f'download and processing of sequences...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    skipped_samples = 0
    for sample_id in tqdm(accession_ids):
        try:
            sample_path = _download_or_get_virus_sample_as_xml(download_sample_dir_path, sample_id)
        except urllib.error.URLError:
            logger.exception(f"Network error while downloading sample {sample_id} can't be solved. This sample 'll be skipped.")
            skipped_samples += 1
            continue
        other_sample = SampleWrapperClass(sample_path, sample_id)
        yield other_sample
    logger.info(f"{skipped_samples} 've been skipped due to network errors.")


# noinspection PyPep8Naming
def _reference_sample_from_organism(samples_query: str, log_with_name: str, SampleWrapperClass=AnyNCBIVNucSample):
    logger.trace(f'importing reference sample of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of refseq...')
    accession_ids = _get_samples_accession_ids(samples_query+' AND srcdb_refseq[Properties]')
    assert len(accession_ids) > 0, f'No reference sequence exists in the sample group identified by the query\n' \
                                   f'<<<' \
                                   f'{samples_query}' \
                                   f'>>>'
    assert len(accession_ids) == 1, f'More than one reference sequence found in the sample group identified by the query\n' \
                                    f'<<<\n' \
                                    f'{samples_query}\n' \
                                    f'>>>\n' \
                                    f'Reference sequence accession ids: {accession_ids}'

    logger.trace(f'download and processing of the reference sample...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    sample_path = _download_or_get_virus_sample_as_xml(download_sample_dir_path, accession_ids[0])
    ref_sample = SampleWrapperClass(sample_path, accession_ids[0])
    return ref_sample


def _download_virus_taxonomy_as_xml(containing_directory: str, taxon_id: int) -> str:
    def count_organisms():
        with Entrez.esearch(db="taxonomy", term=f'{taxon_id}[uid]', rettype='count', retmode="xml") as handle_1:
            response = Entrez.read(handle_1)
        return int(response['Count'])
    how_many_organisms = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, count_organisms)
    assert how_many_organisms == 1, f'The passed taxon id should match exactly one organism in the NCBI Taxonomy DB but\n' \
                                    f'{taxon_id} does not match any.'

    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{os.path.sep}{taxon_id}.xml"
        if not os.path.exists(destination_file_path):
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                    open(destination_file_path, 'w') as f:
                for line in handle:
                    f.write(line)
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


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
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _get_samples_accession_ids(all_samples_query: str) -> List[int]:
    def do():
        # DO PAGINATION
        # total number of sequences
        with Entrez.esearch(db="nuccore",
                            term=f"{all_samples_query}",
                            rettype='count') as handle:
            response = Entrez.read(handle)
            total_records = int(response['Count'])
        # get pages
        accessions_ids = list()
        # noinspection PyPep8Naming
        RECORDS_PER_PAGE = 1000
        page_number = 0
        with tqdm(total=total_records) as progress_bar:
            while total_records > page_number * RECORDS_PER_PAGE:
                with Entrez.esearch(db="nuccore",
                                    term=f"{all_samples_query}",
                                    retmax=RECORDS_PER_PAGE, retstart=page_number * RECORDS_PER_PAGE) as handle:
                    response = Entrez.read(handle)
                for x in response['IdList']:
                    accessions_ids.append(int(x))
                page_number += 1
                progress_bar.update(
                    RECORDS_PER_PAGE if page_number * RECORDS_PER_PAGE < total_records else total_records - (
                                (page_number - 1) * RECORDS_PER_PAGE))
        if len(accessions_ids) != total_records:
            raise IOError('Some of the accession ids were not correctly downloaded')
        return accessions_ids
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _try_n_times(n_times: int, or_wait_secs: int, function: Callable, *args, **kwargs):
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
                return _try_n_times(n_times, or_wait_secs, function, *args, **kwargs)
            else:
                raise e


def main_pipeline_part_3(session: database_tom.Session, sample: AnyNCBIVNucSample, db_sequence_id):
    try:
        # logger.debug(f'callling sequence aligner with args: {db_sequence_id}, <ref_seq>, <seq>, {virus_sequence_chromosome_name}, {annotation_file_path}, {snpeff_db_name}')
        # logger.debug('seq: '+sample.nucleotide_sequence())
        file_path = get_local_folder_for(log_with_name, FileType.Annotations)+str(sample.internal_id())+".pickle"
        if not os.path.exists(file_path):
            annotations_and_nuc_variants = sequence_aligner(
                sample.internal_id(),
                nucleotide_reference_sequence,
                sample.nucleotide_sequence(),
                virus_sequence_chromosome_name,
                annotation_file_path,
                snpeff_db_name)
            with open(file_path, mode='wb') as cache_file:
                pickle.dump(annotations_and_nuc_variants, cache_file, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(file_path, mode='rb') as cache_file:
                annotations_and_nuc_variants = pickle.load(cache_file)
        annotations, nuc_variants = annotations_and_nuc_variants
        for ann in annotations:
            vcm.create_annotation_and_amino_acid_variants(session, db_sequence_id, *ann)
        for nuc in nuc_variants:
            vcm.create_nuc_variants_and_impacts(session, db_sequence_id, nuc)
    except Exception:
        logger.exception(f'exception occurred while working on annotations and nuc_variants of virus sample {sample}. Rollback transaction.')
        raise database_tom.RollbackTransactionWithoutError()


class ImportAnnotationsAndNucVariants(Job):
    def __init__(self, sequence_id, sample: AnyNCBIVNucSample):
        super(ImportAnnotationsAndNucVariants, self).__init__()
        self.sequence_id = sequence_id
        self.sample = sample
        sample.on_before_multiprocessing()

    def execute(self):
        try:
            main_pipeline_part_3(self.session, self.sample, self.sequence_id)
            self.session.commit()
            logger.info(f'pipeline_part_3 of sequence with id {self.sequence_id} imported')
        except database_tom.RollbackTransactionWithoutError:
            database_tom.rollback(self.session)
        except:
            logger.exception(f'unknown exception while running pipeline_part_3 of sequence with id {self.sequence_id}')
            database_tom.rollback(self.session)

    def use_session(self, session):
        self.session = session


class TheWorker(Worker):
    def __init__(self, jobs):
        super(TheWorker, self).__init__(jobs)
        self.worker_session = database_tom.get_session()

    def get_a_job(self) -> Optional[Job]:
        # noinspection PyTypeChecker
        job: ImportAnnotationsAndNucVariants = super().get_a_job()
        if job:
            job.use_session(self.worker_session)
        return job

    def release_resources(self):
        super().release_resources()
        self.worker_session.close()


# noinspection PyPep8Naming
def import_samples_into_vcm_except_annotations_nuc_vars(
        samples_query,
        bind_to_organism_taxon_id: int,
        _log_with_name: str,
        _annotation_file_path: str,
        virus_chromosome_name: str,
        _snpeff_db_name: str,
        SampleWrapperClass=AnyNCBIVNucSample,
        OrganismWrapperClass=AnyNCBITaxon,
        from_sample: Optional[int] = None,
        to_sample: Optional[int] = None):

    global annotation_file_path, virus_sequence_chromosome_name, snpeff_db_name, log_with_name
    annotation_file_path = _annotation_file_path
    virus_sequence_chromosome_name = virus_chromosome_name
    snpeff_db_name = _snpeff_db_name
    log_with_name = _log_with_name
    the_boss = Boss(90, 90, TheWorker)

    def check_user_query():
        """
        Makes sure the user-provided query returns at least one sequence.
        :raise AssertionError if the user query returns 0 sequences.
        """
        # total number of sequences
        with Entrez.esearch(db="nuccore",
                            term=f"{samples_query}",
                            rettype='count') as handle:
            response = Entrez.read(handle)
        total_records = int(response['Count'])
        assert total_records > 0, f'The query provided does not select any sample from NCBI Nucleotide DB. ' \
                                  f'Retry with a different query \n' \
                                  f'<<<\n' \
                                  f'{samples_query}\n' \
                                  f'>>>'
        logger.info(f'The query provided selects {total_records} form NCBI Nucleotide DB.')

    def main_pipeline_part_1():
        """
        Inserts or fetches the organism the user wants to bind the sequences to into the VIRUS table.
        Also it fetches the reference nucleotide sequence for such organism taking it from the DB if already present,
        or extracting it from the samples selected from the user-provided query.
        :return:
        """
        taxonomy_download_dir = get_local_folder_for(log_with_name, FileType.TaxonomyData)
        taxon_file_path = _download_virus_taxonomy_as_xml(taxonomy_download_dir, bind_to_organism_taxon_id)
        organism = OrganismWrapperClass(taxon_file_path)

        def get_virus_id_and_nuc_reference_sequence(session: database_tom.Session):
            # get virus ID
            virus = vcm.get_virus(session, organism)
            if not virus:
                logger.trace(f'inserting new organism {bind_to_organism_taxon_id} into VIRUS table')
                virus = vcm.create_or_get_virus(session, organism)
            # get nuc_reference_sequence
            logger.trace(f'organism txid {bind_to_organism_taxon_id} already in VIRUS table.')
            refseq_row: database_tom.Sequence = vcm.get_reference_sequence_of_virus(session, virus)
            if refseq_row:
                return virus.virus_id, refseq_row.nucleotide_sequence
            else:
                logger.trace('caching the nucleotide sequence of the reference')
                downloaded_ref_sample = _reference_sample_from_organism(samples_query, log_with_name, SampleWrapperClass)
                return virus.virus_id, downloaded_ref_sample.nucleotide_sequence()

        return database_tom.try_py_function(
            get_virus_id_and_nuc_reference_sequence
        )

    def main_pipeline_part_2(session: database_tom.Session, sample: AnyNCBIVNucSample):
        """
        Inserts metadata and returns the SEQUENCE.sequence_id
        :return: SEQUENCE.sequence_id of every sample
        """
        try:
            experiment = vcm.create_or_get_experiment(session, sample)
            host_sample = vcm.create_or_get_host_sample(session, sample)
            sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
            sequence = vcm.create_or_get_sequence(session, sample, virus_id, experiment, host_sample, sequencing_project)
            return sequence.sequence_id
        except Exception as e:
            if str(e).startswith('duplicate key value violates unique constraint "sequence_accession_id_key"'):
                logger.error(f'exception occurred while working on virus sample {sample}: {str(e)}')
            else:
                logger.exception(f'exception occurred while working on virus sample {sample}')
            raise database_tom.RollbackTransactionWithoutError()

    check_user_query()

    global nucleotide_reference_sequence
    virus_id, nucleotide_reference_sequence = main_pipeline_part_1()   # id of the virus in the VCM.Virus table, to which imported sequences will be bound

    logger.info(f'importing the samples identified by the query provided as bound to organism txid{bind_to_organism_taxon_id}')
    database_tom.dispose_db_engine()
    the_boss.wake_up_workers()
    database_tom.re_config_db_engine(False)
    # counter = 30
    try:
        for sample in _all_samples_from_organism(samples_query, log_with_name, SampleWrapperClass, from_sample, to_sample):
        # if sample.length() < 29000: # TODO remove this
        #     continue
        # if counter > 0:             # TODO remove this
        #     counter -= 1
            sequence_id = database_tom.try_py_function(
                main_pipeline_part_2, sample
            )
            if sequence_id:
                # carries out main_pipeline_3 with many processes
                the_boss.assign_job(ImportAnnotationsAndNucVariants(sequence_id, sample))
                # database_tom.try_py_function(
                #     main_pipeline_part_3, sample, sequence_id
                # )
        # else:
        #     break
    finally:
        # pass
        the_boss.stop_workers()



prepared_parameters = {
    'dengue_virus_1': ('txid11053[Organism:exp]', 11053, 'Dengue Virus 1', f'.{sep}annotations{sep}dengue_virus_1.tsv', 'NC_001477', 'dengue_virus_1'),
    'dengue_virus_2': ('txid11060[Organism:exp]', 11060, 'Dengue Virus 2', f'.{sep}annotations{sep}dengue_virus_2.tsv', 'NC_001474', 'dengue_virus_2'),
    'dengue_virus_3': ('txid11069[Organism:exp]', 11069, 'Dengue Virus 3', f'.{sep}annotations{sep}dengue_virus_3.tsv', 'NC_001475', 'dengue_virus_3'),
    'dengue_virus_4': ('txid11070[Organism:exp]', 11070, 'Dengue Virus 4', f'.{sep}annotations{sep}dengue_virus_4.tsv', 'NC_002640', 'dengue_virus_4'),
    'mers': ('txid1335626[Organism:noexp]', 1335626, 'MERS-CoV', f'.{sep}annotations{sep}mers.tsv', 'NC_019843', 'mers'),
    'betacoronavirus_england_1': ('txid1263720[Organism:noexp]', 1263720, 'Betacoronavirus England 1', f'.{sep}annotations{sep}betacoronavirus_england_1.tsv', 'NC_038294', 'betacoronavirus_england_1'),
    'zaire_ebolavirus': ('txid186538[Organism:exp]', 186538, 'Zaire ebolavirus', f'.{sep}annotations{sep}zaire_ebolavirus.tsv', 'NC_002549', 'zaire_ebolavirus'),
    'sudan_ebolavirus': ('txid186540[Organism:exp]', 186540, 'Sudan ebolavirus', f'.{sep}annotations{sep}sudan_ebolavirus.tsv', 'NC_006432', 'sudan_ebolavirus'),
    'reston_ebolavirus': ('txid186539[Organism:exp]', 186539, 'Reston ebolavirus', f'.{sep}annotations{sep}reston_ebolavirus.tsv', 'NC_004161', 'reston_ebolavirus'),
    'bundibugyo_ebolavirus': ('txid565995[Organism:noexp]', 565995, 'Bundibugyo ebolavirus', f'.{sep}annotations{sep}bundibugyo_ebolavirus.tsv', 'NC_014373', 'bundibugyo_ebolavirus'),
    'bombali_ebolavirus': ('txid2010960[Organism:noexp]', 2010960, 'Bombali ebolavirus', f'.{sep}annotations{sep}bombali_ebolavirus.tsv', 'NC_039345', 'bombali_ebolavirus'),
    'tai_forest_ebolavirus': ('txid186541[Organism:exp]', 186541, 'Tai Forest ebolavirus', f'.{sep}annotations{sep}tai_forest_ebolavirus.tsv', 'NC_014372', 'tai_forest_ebolavirus'),
    'new_ncbi_sars_cov_2': ('txid2697049[Organism]', 2697049, 'New NCBI SARS-Cov-2', f'.{sep}annotations{sep}new_ncbi_sars_cov_2.tsv', 'NC_045512', 'new_ncbi_sars_cov_2'),
    'new_ncbi_sars_cov_1': ('txid694009[Organism:noexp] NOT txid2697049[Organism]', 694009, 'New NCBI SARS-Cov-1', f'.{sep}annotations{sep}sars_cov_1.tsv', 'NC_004718', 'sars_cov_1')
}
