from os.path import sep
from collections import Counter, OrderedDict
from data_sources.ncbi_any_virus.host_sample import NCBIHostSample
from pipeline_nuc_variants__annotations__aa import sequence_aligner
import re
import os
import urllib
import stats_module
from queuable_jobs import Job, Boss, Worker
from datetime import datetime
from decimal import Decimal
from typing import Callable, List, Optional, Tuple, Iterator
from data_sources.common_methods_virus import _try_n_times, download_ncbi_taxonomy_as_xml, download_or_get_ncbi_sample_as_xml
from data_sources.common_methods_host_sample import host_taxon_id_from_ncbi_taxon_name, host_taxon_name_from_ncbi_taxon_id
import dateutil.parser as dateparser
import lxml
from tqdm import tqdm
import database_tom
import vcm as vcm
from geo_groups import geo_groups
from xml_helper import text_at_node
from loguru import logger
from lxml import etree
from locations import get_local_folder_for, FileType, remove_file
from Bio import Entrez
import pickle
import cleaning_module
Entrez.email = "Your.Name.Here@example.org"


nucleotide_reference_sequence: Optional[str] = None # initialized elsewhere
annotation_file_path: Optional[str] = None  # initialized elsewhere
virus_sequence_chromosome_name: Optional[str] = None  # initialized elsewhere
snpeff_db_name: Optional[str] = None  # initialized elsewhere
log_with_name: Optional[str] = None
user_provided_reference_query: Optional[str] = None
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
    re_anna_RuL3z = re.compile(r'\D*(\d+)([,|.]?)(\d*).*')
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)
    re_patent_submission_date_lab = re.compile(r'^.* (\d+-\w+-\d+)(?: (.+))?$')
    gender_suggested_by_other_method: Optional[str] = None
    country_suggested_by_other_method: Optional[str] = None
    host_name_suggested_by_other_method:Optional[str] = None
    isolation_source_suggested_by_other_method: Optional[str] = None

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
        self._nuc_seq = None
        self._host_taxon_name = None
        self._host_taxon_id = None
        self._external_host_data = None

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
            # Note: isolate is synonymous of strain and completely different from isolation_source. If it happens
            # to find isolations_source values in this field, then it's an error of the data source.
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
            if strain == 'nasopharyngeal':
                strain = None
                isolation_source_suggested_by_other_method = 'nasopharyngeal'
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

    def coverage(self) -> Optional[int]:
        _input = self._structured_comment(self.sample_xml, 'Coverage')
        if not _input:
            info = None
        else:
            info = self.re_anna_RuL3z.match(_input) # result of regular expression to extract the value may be None
            if not info:
                logger.warning(f'coverage string in sample {self.internal_id()} is {_input} and cannot be parsed. '
                               f'EXPERIMENT_TYPE.coverage set to NULL.\n'
                               f'Full comment content is: {text_at_node(self.sample_xml, ".//INSDSeq_comment", False)}')
            else:
                o1 = info.group(1)  # integer part
                o3 = info.group(3)  # decimal (optional)
                # apply Anna's rules
                output = int(o1)
                if o3 and len(o3) < 3:
                    # then o3 is the decimal part and i moust be rounded to closest integer
                    decimals = round(float(f'0.{o3}'))
                    output += decimals
                return output
        # one last attempt
        if not info and self.external_host_data() is not None:
            return self.external_host_data().coverage()

    def collection_date(self) -> Optional[str]:
        collection_date = text_at_node(self.sample_xml,
                                       '..//INSDQualifier[./INSDQualifier_name/text() = "collection_date"]/INSDQualifier_value',
                                       mandatory=False)
        if collection_date:
            try:
                collection_date = dateparser.parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
            except dateparser._parser.ParserError as e:
                if '/' in collection_date:
                    collection_date = dateparser.parse(collection_date.split('/')[0], default=self.default_datetime).strftime('%Y-%m-%d')
                else:
                    collection_date = None
                logger.warning(f'collection date string"{collection_date}" was parsed as {collection_date}')
        if collection_date is None and self.external_host_data() is not None:
            collection_date = self.external_host_data().collection_date()
        return collection_date

    def submission_date(self) -> Optional[datetime]:
        try:
            return datetime.strptime(self._init_and_get_journal()[0], '%d-%b-%Y')
        except TypeError:
            if self.external_host_data() is not None:
                return self.external_host_data().submission_date()
            else:
                return None

    def isolation_source(self):
        source = text_at_node(self.sample_xml,
                              '..//INSDQualifier[./INSDQualifier_name/text() = "isolation_source"]/INSDQualifier_value',
                              mandatory=False)
        if not source:
            source = self.isolation_source_suggested_by_other_method
        if source:
            source = source.lower()
            # sometimes isolation_source contains completely different data. In that case, the value is
            # useless for the purpose of obtaining the isolation source so we return None
            country = geo_groups.get(source.lower())

            # check if it contains gender info
            if 'female' in source:
                self.gender_suggested_by_other_method = 'female'
                if source == 'female':
                    source = None
            elif 'male' in source:
                self.gender_suggested_by_other_method = 'male'
                if source == 'male':
                    source = None
            elif 'sapiens' in source:
                self.host_name_suggested_by_other_method = 'Homo sapiens'
                source = None
            elif country is not None:
                self.country_suggested_by_other_method = country
                source = None
        # last chance
        if source is None and self.external_host_data() is not None:
            return self.external_host_data().isolation_source()
        else:
            return source

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
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
        if not country and not region and self.external_host_data() is not None:
            country, region, geo_group = self.external_host_data().country__region__geo_group()
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

    def originating_lab(self) -> Optional[str]:
        if self.external_host_data() is not None:
            return self.external_host_data().originating_lab()
        else:
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

    def host_taxon_name(self) -> Optional[str]:
        if not self._host_taxon_name:   # cached host_taxon_name
            if not self.host_name_suggested_by_other_method:
                # then find in sample XML
                host_values = self._init_and_get_host_values()
                host = None
                if host_values:
                    host = host_values[0]
                    try:    # remove stuff in parentheses
                        first_par = host.index('(')  # raises ValueError if ( is not found
                        after_second_par = host.index(')', first_par) + 1
                        output = host[:first_par].strip()
                        if after_second_par < len(host):
                            output = output + ' ' + host[after_second_par:].strip()
                        host = output
                    except ValueError:
                        pass
                # if not foud in sequence XML, look into external host XML file
                if not host and self.external_host_data() is not None:
                    host = self.external_host_data().host_taxon_name()
            else:
                host = self.host_name_suggested_by_other_method
            # correct typos
            if host:
                host = cleaning_module.correct_typos(host)
                self._host_taxon_name = host
            # at this point we may have None or a somehow valid host_taxon_name
            # let's try to uniform this value with its synonyms by searching it in NCBI
            # the NCBI API to search returns only IDs
            txid = host_taxon_id_from_ncbi_taxon_name(self._host_taxon_name)
            self._host_taxon_id = txid  # cache the taxon_id
            ncbi_taxon_name = host_taxon_name_from_ncbi_taxon_id(txid)
            if ncbi_taxon_name:
                self._host_taxon_name = ncbi_taxon_name
        return self._host_taxon_name

    # host taxon id
    def host_taxon_id(self) -> Optional[int]:
        taxon_name = self.host_taxon_name()  # this call may initialize the cached taxon id
        return self._host_taxon_id or host_taxon_id_from_ncbi_taxon_name(taxon_name)


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
            gender = gender or self.gender_suggested_by_other_method      # can be None
        else:
            gender = self.gender_suggested_by_other_method    # can be None
        if not gender and self.external_host_data() is not None:
            gender = self.external_host_data().gender()
        return gender

    def age(self) -> Optional[int]:
        host = self._init_and_get_host_values()
        age = None
        if host:
            age = next(filter(lambda x: 'age' in x, host), None)
            if age:
                age = int(''.join(filter(lambda x: x.isdigit(), age)))
                # age = int(next(filter(lambda x: x.isdigit(), age.split()), None))
        if not age and self.external_host_data() is not None:
            age = self.external_host_data().age()
        return age

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
        self._external_host_data = None

    def has_external_host_data(self):
        nodes = self.sample_xml.xpath('.//INSDSeq_xrefs/INSDXref[./INSDXref_dbname/text() = "BioSample"]')
        if len(nodes) == 1:
            try:
                host_id = text_at_node(nodes[0], './INSDXref_id', True)
                return True
            except AssertionError:
                return False
        else:
            return False

    def external_host_data(self):
        # check cached value first
        if self._external_host_data == -1:      # already initialized but it doesn't have external host data
            return None
        elif self._external_host_data is None:  # never initialized
            # compute value
            nodes = self.sample_xml.xpath('.//INSDSeq_xrefs/INSDXref[./INSDXref_dbname/text() = "BioSample"]')
            if len(nodes) == 1:
                try:
                    host_id = text_at_node(nodes[0], './INSDXref_id', True)
                    self._external_host_data = NCBIHostSample(host_id)
                except AssertionError as e:
                    logger.error('malforned XML file. BioSample attribute without ID')
                    raise e
            else:
                self._external_host_data = -1
                return None
        return self._external_host_data


def _alt_ids_of_selected_sequences(samples_query: str, log_with_name: str) -> List[int]:
    logger.trace(f'importing samples of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of samples...')
    accession_ids = _get_samples_accession_ids(samples_query)

    # if the reference sample is not selected already from samples_query, add it as first sample
    global user_provided_reference_query
    if user_provided_reference_query:
        accession_ids.insert(0, _get_samples_accession_ids(user_provided_reference_query)[0])
    return accession_ids


# noinspection PyPep8Naming
def _download_as_sample_object(alternative_accession_ids, log_with_name: str, SampleWrapperClass=AnyNCBIVNucSample) -> Iterator:
    logger.trace(f'download and processing of sequences...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    network_skipped_samples = 0
    parser_skipped_samples = 0
    other_reasons_skipped_samples = 0
    for sample_id in tqdm(alternative_accession_ids):
        try:
            sample_path = download_or_get_ncbi_sample_as_xml(download_sample_dir_path, sample_id)
        except urllib.error.URLError:
            logger.exception(
                f"Network error while downloading sample {sample_id} can't be solved. This sample 'll be skipped.")
            network_skipped_samples += 1
            #  downloaded file may be incomplete
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            continue
        except KeyboardInterrupt:
            #  downloaded file may be incomplete
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            return
        except TypeError:
            logger.error(f"Error while writing the XML for sample {sample_id}. This sample 'll be skipped.")
            other_reasons_skipped_samples += 1
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            continue
        try:
            other_sample = SampleWrapperClass(sample_path, sample_id)
        except (lxml.etree.XMLSyntaxError, AssertionError):
            logger.exception(f"Error while parsing the sample XML {sample_path}. The file is being deleted, and the sample skipped.")
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            parser_skipped_samples += 1
            continue
        yield other_sample
    logger.info(f"{network_skipped_samples} 've been skipped due to network errors.\n"
                f"{parser_skipped_samples} 've been skipped due to XML parser errors.\n"
                f"{other_reasons_skipped_samples} 've been skipped because of other errors.")


# noinspection PyPep8Naming
def _reference_sample_from_organism(samples_query: str, log_with_name: str, SampleWrapperClass=AnyNCBIVNucSample):
    logger.trace(f'importing reference sample of organism "{log_with_name}"')

    logger.trace(f'getting accession ids of refseq...')
    reference_query = samples_query+' AND srcdb_refseq[Properties]'
    accession_ids = _get_samples_accession_ids(reference_query)

    def test_query():
        global user_provided_reference_query
        nonlocal accession_ids, reference_query
        if len(accession_ids) == 0:
            logger.warning(f'No reference sequence exists in the sample group identified by the query\n'
                           f'<<<\n'
                           f'{reference_query}\n'
                           f'>>>\n'
                           f'Type a valid query selecting only one sample to be used as reference sample or press Ctrl+C to abort')
            user_provided_reference_query = input().strip()
            reference_query = user_provided_reference_query
            accession_ids = _get_samples_accession_ids(reference_query)
            test_query()
        if len(accession_ids) != 1:
            logger.warning(f'More than one reference sequence found in the sample group identified by the query\n'
                           f'<<<\n'
                           f'{reference_query}\n'
                           f'>>>\n'
                           f'Reference sequence accession ids: {accession_ids}\n'
                           f'Type a valid query selecting only one sample to be used as reference sample or press Ctrl+C to abort')
            user_provided_reference_query = input().strip()
            reference_query = user_provided_reference_query
            accession_ids = _get_samples_accession_ids(reference_query)
            test_query()

    test_query()
    logger.trace(f'download and processing of the reference sample...')
    download_sample_dir_path = get_local_folder_for(source_name=log_with_name, _type=FileType.SequenceOrSampleData)
    try:
        sample_path = download_or_get_ncbi_sample_as_xml(download_sample_dir_path, accession_ids[0])
    except Exception as e:
        logger.error('Error while downloading or getting the reference sample XML. Partial file will be deleted')
        remove_file(f"{download_sample_dir_path}{accession_ids[0]}.xml")
        raise e
    ref_sample = SampleWrapperClass(sample_path, accession_ids[0])
    return ref_sample


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
        RECORDS_PER_PAGE = 5000
        page_number = 0
        import time
        with tqdm(total=total_records) as progress_bar:
            while total_records > page_number * RECORDS_PER_PAGE:
                with Entrez.esearch(db="nuccore",
                                    term=f"{all_samples_query}",
                                    retmax=RECORDS_PER_PAGE, retstart=page_number * RECORDS_PER_PAGE) as handle:
                    response = Entrez.read(handle)
                for x in response['IdList']:
                    accessions_ids.append(int(x))
                page_number += 1
                time.sleep(2)
                progress_bar.update(
                    RECORDS_PER_PAGE if page_number * RECORDS_PER_PAGE < total_records else total_records - (
                                (page_number - 1) * RECORDS_PER_PAGE))
        if len(accessions_ids) != total_records:
            raise IOError('Some of the accession ids were not correctly downloaded')
        return accessions_ids
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _deltas(virus_id, id_remote_samples):
    """
    :param virus_id: virus_id from the DB of the virus being imported
    :param id_remote_samples: alternative accession id of the sequences avbailable from remote
    :return: three sets containing the alternative accession ids of
        1. sequences that are unchanged (local == remote)
        2. outdated sequences (to be removed)
        3. sequences new and to be imported
    """
    id_local_samples = set([int(i) for i in database_tom.try_py_function(
        vcm.sequence_alternative_accession_ids, virus_id, ['GenBank', 'RefSeq']
    )])
    id_outdated_sequences = id_local_samples - id_remote_samples
    id_new_sequences = id_remote_samples - id_local_samples
    id_unchanged_sequences = id_local_samples - id_outdated_sequences
    logger.info(f'\n'
                f'# Sequences from remote source: {len(id_remote_samples)}. Of which\n'
                f'# {len(id_unchanged_sequences)} present also locally and unchanged\n'
                f'# {len(id_new_sequences)} never seen before.\n'
                f'# Sequences from local source: {len(id_local_samples)}. Of which\n'
                f'# {len(id_outdated_sequences)} outdated and must be removed from local')
    return id_unchanged_sequences, id_outdated_sequences, id_new_sequences


def main_pipeline_part_3(session: database_tom.Session, sample: AnyNCBIVNucSample, db_sequence_id):
    file_path = get_local_folder_for(log_with_name, FileType.Annotations) + str(sample.internal_id()) + ".pickle"
    try:
        # logger.debug(f'callling sequence aligner with args: {db_sequence_id}, <ref_seq>, <seq>, {virus_sequence_chromosome_name}, {annotation_file_path}, {snpeff_db_name}')
        # logger.debug('seq: '+sample.nucleotide_sequence())
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
    except KeyboardInterrupt:
        # cached file may be incomplete
        remove_file(file_path)
        raise database_tom.Rollback()
    except Exception:
        logger.exception(f'exception occurred while working on annotations and nuc_variants of virus sample '
                         f'{sample.internal_id()}. Rollback transaction.')
        raise database_tom.Rollback()


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
            stats_module.completed_sample(self.sample.alternative_accession_number())
        except database_tom.Rollback:
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
def import_samples_into_vcm(
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
    the_boss = Boss(70, 70, TheWorker)

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
        taxon_file_path = download_ncbi_taxonomy_as_xml(taxonomy_download_dir, bind_to_organism_taxon_id)
        organism = OrganismWrapperClass(taxon_file_path)

        def get_virus_id_and_nuc_reference_sequence(session: database_tom.Session):
            # get virus ID
            virus = vcm.get_virus(session, organism)
            if virus:
                _virus_id = virus.virus_id
            else:
                logger.trace(f'inserting new organism {bind_to_organism_taxon_id} into VIRUS table')
                _virus_id = vcm.create_or_get_virus(session, organism)
            # get nuc_reference_sequence
            logger.trace(f'organism txid {bind_to_organism_taxon_id} already in VIRUS table.')
            refseq_row: database_tom.Sequence = vcm.get_reference_sequence_of_virus(session, _virus_id)
            if refseq_row:
                return _virus_id, refseq_row.nucleotide_sequence
            else:
                logger.trace('caching the nucleotide sequence of the reference')
                downloaded_ref_sample = _reference_sample_from_organism(samples_query, log_with_name, SampleWrapperClass)
                return _virus_id, downloaded_ref_sample.nucleotide_sequence()

        return database_tom.try_py_function(
            get_virus_id_and_nuc_reference_sequence
        )

    def main_pipeline_part_2(session: database_tom.Session, sample: AnyNCBIVNucSample):
        """
        Inserts metadata and returns the SEQUENCE.sequence_id
        :return: SEQUENCE.sequence_id of every sample
        """
        try:
            experiment_id = vcm.create_or_get_experiment(session, sample)
            host_specie_id = vcm.create_or_get_host_specie(session, sample)
            host_sample_id = vcm.create_or_get_host_sample(session, sample, host_specie_id)
            sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
            sequence = vcm.create_and_get_sequence(session, sample, virus_id, experiment_id, host_sample_id, sequencing_project_id)
            return sequence.sequence_id
        except ValueError as e:
            if str(e).endswith('missing nucleotide seq.'):
                logger.warning(f'Sample {sample.internal_id()} is given without nucleotide sequence. Sample import skipped. File deleted')
                remove_file(f"{get_local_folder_for(log_with_name, FileType.SequenceOrSampleData)}{sample.alternative_accession_number()}.xml")
            else:
                logger.exception(f"exception occurred while working on virus sample {sample.internal_id()}. Sample won't be imported")
            raise database_tom.Rollback()
        except AssertionError:
            logger.exception(f"Sample {sample.internal_id()} may have a corrupted XML. The XML will be deleted. Sample won't be imported")
            remove_file(f"{get_local_folder_for(log_with_name, FileType.SequenceOrSampleData)}{sample.alternative_accession_number()}.xml")
            raise database_tom.Rollback()
        except KeyboardInterrupt as e:
            logger.info(f"import of sample {sample.internal_id()} interrupted. Sample won't be imported")
            raise e  # i.e. it is handled outside but I don't want it to be logged as an unexpected error
        except Exception as e:
            logger.exception(f"exception occurred while working on virus sample {sample.internal_id()}. Sample won't be imported")
            raise database_tom.Rollback()

    check_user_query()

    global nucleotide_reference_sequence

    # handle special case of SARS-Cov-1: The reference sample is not included in the sample set.
    if log_with_name == 'New NCBI SARS-Cov-1':
        previous_samples_query = samples_query
        samples_query = 'txid694009[Organism] NOT txid2697049[Organism]'
        virus_id, nucleotide_reference_sequence = main_pipeline_part_1()
        samples_query = previous_samples_query
    else:
        virus_id, nucleotide_reference_sequence = main_pipeline_part_1()  # id of the virus in the VCM.Virus table, to which imported sequences will be bound

    logger.info(f'importing the samples identified by the query provided as bound to organism txid{bind_to_organism_taxon_id}')

    # update last import date
    database_tom.try_py_function(vcm.update_db_metadata, virus_id, 'GenBank')

    # find outdated and new samples from source
    id_all_current_sequences = set(_alt_ids_of_selected_sequences(samples_query, log_with_name))
    _, id_outdated_sequences, id_new_sequences = _deltas(virus_id, id_all_current_sequences)

    # remove outdated sequences
    for alt_seq_acc_id in id_outdated_sequences:
        database_tom.try_py_function(vcm.remove_sequence_and_meta, None, str(alt_seq_acc_id))

    # select range
    id_new_sequences = sorted(list(id_new_sequences))
    if from_sample is not None and to_sample is not None:
        id_new_sequences = id_new_sequences[from_sample:to_sample]

    stats_module.schedule_samples(
        stats_module.StatsBasedOnIds([str(x) for x in id_new_sequences], False, virus_id, ['GenBank', 'RefSeq']))

    # prepare multiprocessing
    database_tom.dispose_db_engine()
    the_boss.wake_up_workers()
    database_tom.re_config_db_engine(False)

    logger.info('begin importing new sequences')
    try:
        for sample in _download_as_sample_object(id_new_sequences, log_with_name, SampleWrapperClass):
            sequence_id = database_tom.try_py_function(
                main_pipeline_part_2, sample
            )
            if sequence_id:
                # carries out main_pipeline_3 with many processes
                the_boss.assign_job(ImportAnnotationsAndNucVariants(sequence_id, sample))
                logger.debug(f'new job queued. waiting jobs: {the_boss.number_of_waiting_jobs()}')
                # database_tom.try_py_function(
                #     main_pipeline_part_3, sample, sequence_id
                # )
    except KeyboardInterrupt:
        # When using multiprocessing, each process receives KeyboardInterrupt. Child process should take care of clean up.
        # Child process should not raise themselves KeyboardInterrupt
        the_boss.discard_left_jobs()
        logger.info('Execution aborted by the user. Cancelling waiting tasks...')
    except Exception as e:
        logger.error("AN EXCEPTION CAUSED THE LOOP TO TERMINATE. Workers will be terminated.")
        raise e
    finally:
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
# QUERY TO SELECT THE REFERENCE OF SARS COV 1
# txid694009[Organism] NOT txid2697049[Organism] AND srcdb_refseq[Properties]
