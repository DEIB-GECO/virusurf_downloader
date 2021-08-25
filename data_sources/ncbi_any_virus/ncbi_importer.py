import sys
from collections import Counter, OrderedDict
from time import sleep
from data_sources.ncbi_any_virus.host_sample import NCBIHostSample
from nuc_aa_pipeline import sequence_aligner
import re
import os
import urllib
import stats_module
from queuable_tasks import Task, TaskManager, Worker
from datetime import datetime
from decimal import Decimal
from typing import Optional, Tuple, Iterator
from data_sources.ncbi_services import host_taxon_id_from_ncbi_taxon_name, host_taxon_name_from_ncbi_taxon_id, \
    download_ncbi_taxonomy_as_xml, download_or_get_ncbi_sample_as_xml, get_samples_accession_ids, read_config_params
import dateutil.parser as dateparser
import lxml
from tqdm import tqdm
from vcm import vcm as vcm
from geo_groups import geo_groups
from xml_helper import text_at_node
from loguru import logger
from lxml import etree
from locations import get_local_folder_for, FileType, remove_file
from Bio import Entrez
import pickle
import data_cleaning_module
from db_config import read_db_import_configuration as db_import_config, database
from data_sources.ncbi_any_virus import settings
from logger_settings import send_message
from http.client import HTTPException
from multiprocessing import Lock
import signal


DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = 30
snpeff_semaphore = Lock()


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
    host_name_suggested_by_other_method: Optional[str] = None
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
        """Used when calling print or similar on this class"""
        return f'AnyNCBIVNucSample:{self.internal_id()}.xml'

    def __repr__(self):
        """Used in stacktraces"""
        return f'AnyNCBIVNucSample:{self.internal_id()}.xml'

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
                isolates = [text_at_node(x, '.', mandatory=True) for x in isolates]
                fallback = isolates[0]
                logger.warning(f'Isolate strings are {isolates}. Only {fallback} was retained as strain')
                return fallback
            if strain == 'nasopharyngeal':
                strain = None
                self.isolation_source_suggested_by_other_method = 'nasopharyngeal'
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
            reference_length = settings.length_of_nuc_reference_sequence
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
            info = self.re_anna_RuL3z.match(_input)  # result of regular expression to extract the value may be None
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

    def collection_date(self) -> Tuple[Optional[str], Optional[int]]:
        collection_date = text_at_node(self.sample_xml,
                                       '..//INSDQualifier[./INSDQualifier_name/text() = "collection_date"]/INSDQualifier_value',
                                       mandatory=False)
        if collection_date:
            if '/' in collection_date:
                collection_date = collection_date.split('/')[0]
            try:
                collection_date = collection_date.strip()
                well_formatted_coll_date = dateparser.parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
                # find precision (year = 0, month = 1, day = 2
                if '-' in collection_date:
                    precision = collection_date.count('-')
                elif len(collection_date) == 4:
                    precision = 0
                elif '/' in collection_date:
                    precision = collection_date.count('/')
                elif '\\' in collection_date:
                    precision = collection_date.count('\\')
                elif ' ' in collection_date:
                    precision = collection_date.count(' ')
                else:
                    raise AssertionError(
                        f'Unable to parse date string {collection_date}. Unexpected separator in date string or no '
                        f'separator.')
            except dateparser._parser.ParserError as e:
                if self.external_host_data() is not None:
                    well_formatted_coll_date, precision = self.external_host_data().collection_date()
                else:
                    logger.error(f'Unable to parse date from string {collection_date}. Format not recognised.')
                    raise e
        elif self.external_host_data() is not None:
            well_formatted_coll_date, precision = self.external_host_data().collection_date()
        else:
            well_formatted_coll_date, precision = None, None
        return well_formatted_coll_date, precision

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

    def province__region__country__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        # assign country and region as indicated into the "country" XML node or as suggested by other class methods
        node = text_at_node(self.sample_xml,
                            '..//INSDQualifier[./INSDQualifier_name/text() = "country"]/INSDQualifier_value',
                            mandatory=False)
        if node:
            node = node.split(":")
            country = node[0] or self.country_suggested_by_other_method     # can be None
            region = node[1].strip() if len(node) > 1 else None
        else:
            country = self.country_suggested_by_other_method    # can be None
            region = None
        province = None

        # find geo_group based on country
        geo_group = geo_groups.get(country.lower()) if country else None

        # (2nd attempt:) look for externally referenced data
        if not country and not region and self.external_host_data() is not None:
            province, region, country, geo_group = self.external_host_data().province__region__country__geo_group()

        # apply special corrections for USA counties
        if country and country.strip().upper() == 'USA':
            region = data_cleaning_module.correct_usa_regions(region)

        return province, region, country, geo_group

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
                host = data_cleaning_module.correct_typos(host)
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
            age = next(filter(lambda x: 'age' in x and 'sewage' not in x.lower(), host), None)
            if age:
                age = int(''.join(filter(lambda x: x.isdigit(), age)))
                # age = int(next(filter(lambda x: x.isdigit(), age.split()), None))
        if not age and self.external_host_data() is not None:
            age = self.external_host_data().age()
        return age

    def database_source(self) -> str:
        return 'RefSeq' if self.is_reference() else 'GenBank'

    def bioproject_id(self):
        bioproject_id = self.sample_xml.xpath('.//INSDXref[./INSDXref_dbname/text() = "BioProject"]')
        # xpath returns a list also for single nodes
        if len(bioproject_id) > 0:
            first_id = text_at_node(bioproject_id[0], './INSDXref_id', False)
            if len(bioproject_id) > 1:
                # check if they are equal
                for other_node in bioproject_id[1:]:
                    other_id = text_at_node(other_node, './INSDXref_id', False)
                    if other_id != first_id:
                        raise AssertionError(f'Multiple distinct bioproject ID found at path '
                                             f'.//INSDXref[./INSDXref_dbname/text() = "BioProject"]')
            return first_id
        else:
            return None

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
                    host_data_download_dir = get_local_folder_for(settings.generated_dir_name, FileType.HostData)
                    self._external_host_data = NCBIHostSample(host_id, host_data_download_dir)
                except AssertionError as e:
                    logger.error('malforned XML file. BioSample attribute without ID')
                    raise e
            else:
                self._external_host_data = -1
                return None
        return self._external_host_data

    def remove_external_host_data_file(self):
        if self._external_host_data is not None and self._external_host_data != -1:
            remove_file(self._external_host_data.file_path)


# noinspection PyPep8Naming
def _download_as_sample_object(alternative_accession_ids, SampleWrapperClass=AnyNCBIVNucSample) -> Iterator:
    logger.info(f'download and processing of sequences...')
    download_sample_dir_path = get_local_folder_for(settings.generated_dir_name, _type=FileType.SequenceOrSampleData)
    network_skipped_samples = 0
    parser_skipped_samples = 0
    other_reasons_skipped_samples = 0
    for sample_id in tqdm(alternative_accession_ids):
        try:
            sample_path = download_or_get_ncbi_sample_as_xml(download_sample_dir_path, sample_id)
            other_sample = SampleWrapperClass(sample_path, sample_id)
        except (urllib.error.URLError, HTTPException):
            logger.exception(
                f"Network error while downloading sample {sample_id} can't be solved. This sample 'll be skipped.")
            network_skipped_samples += 1
            #  downloaded file may be incomplete
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            continue
        except KeyboardInterrupt:
            #  downloaded file may be incomplete
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            raise
        except TypeError:
            logger.error(f"Error while writing the XML for sample {sample_id}. This sample 'll be skipped.")
            other_reasons_skipped_samples += 1
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            continue
        except (lxml.etree.XMLSyntaxError, AssertionError):
            logger.exception(f"Error while parsing the sample XML {sample_path}. The file is being deleted, and the sample skipped.")
            remove_file(file_path=f"{download_sample_dir_path}{os.path.sep}{sample_id}.xml")
            parser_skipped_samples += 1
            continue
        yield other_sample
    logger.info(f"\n\t{network_skipped_samples} 've been skipped due to network errors.\n"
                f"\t{parser_skipped_samples} 've been skipped due to XML parser errors.\n"
                f"\t{other_reasons_skipped_samples} 've been skipped because of other errors.")


# noinspection PyPep8Naming
def _reference_sample_from_organism(SampleWrapperClass=AnyNCBIVNucSample):
    logger.trace(f'importing reference sample of organism "{settings.log_with_name}"')

    logger.trace(f'getting accession ids of refseq...')
    reference_query = settings.reference_sample_query    # initializes a temporary variable
    user_provided_reference_query: Optional[str] = None  # stores user input (if necessary) between executions of test_query
    accession_ids = get_samples_accession_ids(settings.reference_sample_query)

    def check_query():
        nonlocal accession_ids, reference_query, user_provided_reference_query
        if len(accession_ids) == 0:
            logger.warning(f'No reference sequence exists in the sample group identified by the query\n'
                           f'<<<\n'
                           f'{reference_query}\n'
                           f'>>>\n'
                           f'Type a valid query selecting only one sample to be used as reference sample or press Ctrl+C to abort')
            user_provided_reference_query = input().strip()
            reference_query = user_provided_reference_query
            accession_ids = get_samples_accession_ids(reference_query)
            check_query()
        if len(accession_ids) != 1:
            logger.warning(f'More than one reference sequence found in the sample group identified by the query\n'
                           f'<<<\n'
                           f'{reference_query}\n'
                           f'>>>\n'
                           f'Reference sequence accession ids: {accession_ids}\n'
                           f'Type a valid query selecting only one sample to be used as reference sample or press Ctrl+C to abort')
            user_provided_reference_query = input().strip()
            reference_query = user_provided_reference_query
            accession_ids = get_samples_accession_ids(reference_query)
            check_query()

    check_query()    # if reference_sample_query doesn't provide one sample, it asks user input until a valid query is given
    logger.trace(f'download and processing of the reference sample...')
    download_sample_dir_path = get_local_folder_for(source_name=settings.generated_dir_name, _type=FileType.SequenceOrSampleData)
    try:
        sample_path = download_or_get_ncbi_sample_as_xml(download_sample_dir_path, accession_ids[0])
    except Exception as e:
        logger.error('Error while downloading or getting the reference sample XML. Partial file will be deleted')
        remove_file(f"{download_sample_dir_path}{accession_ids[0]}.xml")
        raise e
    ref_sample = SampleWrapperClass(sample_path, accession_ids[0])
    return ref_sample


def _deltas(virus_id, id_remote_samples):
    """
    :param virus_id: virus_id from the DB of the virus being imported
    :param id_remote_samples: alternative accession id of the sequences avbailable from remote
    :return: three sets containing the alternative accession ids of
        1. sequences that are unchanged (local == remote)
        2. outdated sequences (to be removed)
        3. sequences new and to be imported
    """
    id_local_samples = set([int(i) for i in database.try_py_function(
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
                f'# {len(id_outdated_sequences)} outdated and must be removed from local.\n'
                f'Check deltas. The program will resume in 20 sec.')
    try:
        sleep(20)
    except KeyboardInterrupt:
        sys.exit(0)
    return id_unchanged_sequences, id_outdated_sequences, id_new_sequences


def main_pipeline_part_3(session: database.Session, sample: AnyNCBIVNucSample, db_sequence_id, semaphore: Lock):
    file_path = get_local_folder_for(settings.generated_dir_name, FileType.Annotations) + str(sample.internal_id()) + ".pickle"
    try:
        # logger.debug(f'callling sequence aligner with args: {db_sequence_id}, <ref_seq>, <seq>, {virus_sequence_chromosome_name}, {annotation_file_path}, {snpeff_db_name}')
        # logger.debug('seq: '+sample.nucleotide_sequence())
        if not os.path.exists(file_path):
            annotations_and_nuc_variants = sequence_aligner(
                sample.internal_id(),
                settings.nucleotide_reference_sequence,
                sample.nucleotide_sequence(),
                settings.chromosome_name,
                settings.annotation_file_path,
                settings.snpeff_db_name,
                semaphore)
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
        raise database.Rollback()
    except Exception:
        logger.exception(f'exception occurred while working on annotations and nuc_variants of virus sample '
                         f'{sample.internal_id()}. Rollback transaction.')
        remove_file(file_path)
        raise database.Rollback()


class ImportAnnotationsAndNucVariants(Task):
    def __init__(self, sequence_id, sample: AnyNCBIVNucSample):
        super(ImportAnnotationsAndNucVariants, self).__init__()
        self.sequence_id = sequence_id
        self.sample = sample
        sample.on_before_multiprocessing()

    def execute(self):
        try:
            main_pipeline_part_3(self.session, self.sample, self.sequence_id, self.snpeff_semaphore)
            self.session.commit()
            stats_module.completed_sample(self.sample.alternative_accession_number())
        except database.Rollback:
            database.rollback(self.session)
        except:
            logger.exception(f'unknown exception while running pipeline_part_3 of sequence with id {self.sequence_id}')
            database.rollback(self.session)

    def use_session(self, session):
        self.session = session

    def use_semaphore(self, semaphore):
        self.snpeff_semaphore = semaphore


class TheWorker(Worker):
    def __init__(self, tasks):
        super(TheWorker, self).__init__(tasks)
        self.worker_session = database.get_session()
        self.semaphore = snpeff_semaphore

    def pick_up_task(self) -> Optional[Task]:
        # noinspection PyTypeChecker
        task: ImportAnnotationsAndNucVariants = super().pick_up_task()
        if task:
            task.use_session(self.worker_session)
            task.use_semaphore(self.semaphore)
        return task

    def release_resources(self):
        super().release_resources()
        self.worker_session.close()


# noinspection PyPep8Naming
def import_samples_into_vcm(source_name: str, SampleWrapperClass=AnyNCBIVNucSample, OrganismWrapperClass=AnyNCBITaxon,
                            from_sample: Optional[int] = None, to_sample: Optional[int] = None):
    interrupted_by_user = False
    # initialize globals
    settings.initialize_settings(*tuple(settings.known_settings[source_name].values()))
    # setup db connection
    db_params: dict = db_import_config.get_database_config_params()
    database.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])
    # prepare multiprocessor
    task_manager = TaskManager(30, 22, TheWorker)

    def check_user_query():
        """
        Makes sure the user-provided query returns at least one sequence.
        :raise AssertionError if the user query returns 0 sequences.
        """
        entrez_config = read_config_params()
        # total number of sequences
        with Entrez.esearch(db="nuccore",
                            term=f"({settings.non_reference_samples_query}) OR ({settings.reference_sample_query})",
                            rettype='count', tool=entrez_config[0], email=entrez_config[1], api_key=entrez_config[2]) as handle:
            response = Entrez.read(handle)
        total_records = int(response['Count'])
        assert total_records > 0, f'The query provided does not select any sample from NCBI Nucleotide DB. ' \
                                  f'Retry with a different query \n' \
                                  f'<<<\n' \
                                  f'{settings.non_reference_samples_query}\n' \
                                  f'>>>'
        logger.info(f'The query provided selects {total_records} form NCBI Nucleotide DB.')

    def main_pipeline_part_1():
        """
        Inserts or fetches the organism the user wants to bind the sequences to into the VIRUS table.
        Also it fetches the reference nucleotide sequence for such organism taking it from the DB if already present,
        or taking it from the user-provided query.
        :return:
        """
        taxonomy_download_dir = get_local_folder_for(settings.generated_dir_name, FileType.TaxonomyData)
        taxon_file_path = download_ncbi_taxonomy_as_xml(taxonomy_download_dir,  settings.virus_taxon_id)
        organism = OrganismWrapperClass(taxon_file_path)

        def get_virus_id_and_nuc_reference_sequence(session: database.Session):
            # get virus ID
            virus = vcm.get_virus(session, organism)
            if virus:
                _virus_id = virus.virus_id
            else:
                logger.trace(f'inserting new organism { settings.virus_taxon_id} into VIRUS table')
                _virus_id = vcm.create_or_get_virus(session, organism)
            # get nuc_reference_sequence
            logger.trace(f'organism txid { settings.virus_taxon_id} already in VIRUS table.')
            refseq_objects = vcm.get_reference_sequence_of_virus(session, _virus_id)
            if refseq_objects:
                ref_seq_obj, ref_nuc_seq_obj = refseq_objects
                return _virus_id, ref_nuc_seq_obj.nucleotide_sequence, ref_seq_obj.accession_id
            else:
                logger.trace('caching the nucleotide sequence of the reference')
                downloaded_ref_sample = _reference_sample_from_organism(SampleWrapperClass)
                return _virus_id, downloaded_ref_sample.nucleotide_sequence(), downloaded_ref_sample.primary_accession_number()

        return database.try_py_function(get_virus_id_and_nuc_reference_sequence)

    # noinspection PyTypeChecker
    def main_pipeline_part_2(session: database.Session, _sample: AnyNCBIVNucSample):
        """
        Inserts metadata and returns the SEQUENCE.sequence_id
        :return: SEQUENCE.sequence_id of every sample
        """
        try:
            experiment_id = vcm.create_or_get_experiment(session, _sample)
            host_specie_id = vcm.create_or_get_host_specie(session, _sample)
            host_sample_id = vcm.create_or_get_host_sample(session, _sample, host_specie_id)
            sequencing_project_id = vcm.create_or_get_sequencing_project(session, _sample)
            sequence, nucleotide_seq = vcm.create_and_get_sequence(session, _sample, virus_id, experiment_id, host_sample_id, sequencing_project_id)
            vcm.DBCache.commit_changes()
            return sequence.sequence_id
        except ValueError as e:
            if str(e).endswith('missing nucleotide seq.'):
                logger.warning(f'Sample {_sample.internal_id()} is given without nucleotide sequence. Sample import skipped. Sample and host files will be deleted.')
            elif str(e).startswith(f"can't find the biosample numeric id for biosample"):
                logger.warning(f'Sample {_sample.internal_id()} is given with invalid external host information. Sample import skipped. Sample and host files will be deleted.')
            else:
                logger.exception(f"exception occurred while working on virus sample {_sample.internal_id()}. Sample won't be imported. Sample import skipped. Sample and host files will be deleted.")
            _sample.remove_external_host_data_file()
            remove_file(f"{get_local_folder_for(settings.generated_dir_name, FileType.SequenceOrSampleData)}{_sample.alternative_accession_number()}.xml")
            vcm.DBCache.rollback_changes()
            raise database.Rollback()
        except AssertionError as e:
            if str(e).endswith('distinct bioproject ID found at path .//INSDXref[./INSDXref_dbname/text() = "BioProject"]'):
                logger.warning(f"Sample {_sample.internal_id()} has >1 bioprojects ID. Sample won't be imported. Sample and host files will be deleted.")
            else:
                logger.exception(f"Sample {_sample.internal_id()} may have a corrupted XML. Sample won't be imported. Sample and host files will be deleted.")
            _sample.remove_external_host_data_file()
            remove_file(f"{get_local_folder_for(settings.generated_dir_name, FileType.SequenceOrSampleData)}{_sample.alternative_accession_number()}.xml")
            vcm.DBCache.rollback_changes()
            raise database.Rollback()
        except KeyboardInterrupt as e:
            logger.info(f"import of sample {_sample.internal_id()} interrupted. Sample won't be imported")
            vcm.DBCache.rollback_changes()
            database.rollback(session)
            raise e  # i.e. it is handled outside but I don't want it to be logged as an unexpected error
        except Exception as e:
            logger.exception(f"exception occurred while working on virus sample {_sample.internal_id()}. Sample won't be imported. Sample and host files will be deleted.")
            _sample.remove_external_host_data_file()
            remove_file(f"{get_local_folder_for(settings.generated_dir_name, FileType.SequenceOrSampleData)}{_sample.alternative_accession_number()}.xml")
            vcm.DBCache.rollback_changes()
            raise database.Rollback()

    check_user_query()

    # id of the virus in the VCM.Virus table, to which imported sequences will be bound
    virus_id, nucleotide_reference_sequence, ref_acc_id = main_pipeline_part_1()
    settings.set_nucleotide_refseq(nucleotide_reference_sequence)    # set as globally available var

    logger.info(f'importing the samples identified by the query provided as bound to organism txid{ settings.virus_taxon_id}')

    # update last import date
    database.try_py_function(vcm.update_db_metadata, virus_id, 'GenBank')

    # find outdated and new samples from source
    id_all_current_sequences = set(get_samples_accession_ids(f"({settings.non_reference_samples_query}) OR ({settings.reference_sample_query})"))
    _, id_outdated_sequences, id_new_sequences = _deltas(virus_id, id_all_current_sequences)

    # select range
    id_new_sequences = sorted(list(id_new_sequences))
    if from_sample is not None and to_sample is not None:
        id_new_sequences = id_new_sequences[from_sample:to_sample]

    # create pipeline_event (will be inserted later)
    pipeline_event = database.PipelineEvent(
        event_date=datetime.now().strftime("%Y-%m-%d"),
        event_name=f'GenBank {source_name} sequences update',
        removed_items=len(id_outdated_sequences),
        changed_items=0,
        added_items=len(id_new_sequences),  # may eventually change if some sequence are not imported
    )

    stats_module.schedule_samples(
        stats_module.StatsBasedOnIds([str(x) for x in id_new_sequences], False, virus_id, ['GenBank', 'RefSeq']))

    # remove outdated sequences
    id_outdated_sequences = [str(x) for x in id_outdated_sequences]
    database.try_py_function(vcm.remove_sequence_and_meta_list, alternative_sequence_accession_id=id_outdated_sequences)
    stats_module.removed_samples(id_outdated_sequences)

    # prepare multiprocessing
    database.dispose_db_engine()
    logger.warning("Worker processes will ignore SIGINT. The parent process must handle this condition.")
    # temporary disables SIGINT handler in parent and restore it after the fork
    default_ctrl_c_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    task_manager.wake_up_workers()
    signal.signal(signal.SIGINT, default_ctrl_c_handler)
    database.re_config_db_engine(False)

    logger.info('begin importing new sequences')
    vcm.DBCache.commit_changes()
    try:
        for sample in _download_as_sample_object(id_new_sequences, SampleWrapperClass):
            sequence_id = database.try_py_function(
                main_pipeline_part_2, sample
            )
            if sequence_id:
                # carries out main_pipeline_3 with many processes
                task_manager.assign_task(ImportAnnotationsAndNucVariants(sequence_id, sample))
                # logger.debug(f'new job queued. waiting jobs: {task_manager.number_of_waiting_tasks()}')
                # ... or do it one by one
                # database.try_py_function(
                #     main_pipeline_part_3, sample, sequence_id
                # )
    except KeyboardInterrupt:
        # When using multiprocessing, each process receives KeyboardInterrupt. Children processes should take care of clean up.
        # Children processes should not raise themselves KeyboardInterrupt
        interrupted_by_user = True
        logger.info('Execution aborted by the user. Cancelling waiting tasks...')
        task_manager.discard_waiting_tasks()
    except Exception as e:
        logger.error("AN EXCEPTION CAUSED THE LOOP TO TERMINATE. Workers will be terminated.")
        task_manager.discard_waiting_tasks()
        raise e
    finally:
        # ignore Ctrl+C in future
        logger.warning("SIGINT is being ignored during DB finalization process.")
        signal.signal(signal.SIGINT, signal.SIG_IGN)

        task_manager.stop_workers()

        # remove leftovers of failed samples
        try:
            metadata_samples_to_remove: set = stats_module.get_scheduled_not_completed()
            if len(metadata_samples_to_remove) > 100 and not interrupted_by_user:
                send_message(f"NCBI importer can have a bug. {len(metadata_samples_to_remove)} out of "
                             f"{len(id_new_sequences)} failed.")
            pipeline_event.added_items = pipeline_event.added_items - len(metadata_samples_to_remove)
            if len(metadata_samples_to_remove) > 0:
                logger.info(f"Removing metadata leftovers of imports that failed during variant/annotation calling or metadata"
                            f" ({len(metadata_samples_to_remove)} samples)")

                metadata_samples_to_remove_as_string: list = [str(x) for x in metadata_samples_to_remove]
                logger.trace("Alternative accession id of failed imports:\n"
                             f"{metadata_samples_to_remove_as_string}")
                logger.info("Deleting leftovers in database")
                database.try_py_function(vcm.remove_sequence_and_meta_list,
                                         alternative_sequence_accession_id=metadata_samples_to_remove_as_string)
                logger.info("Deleting XML sequence files")
                for x in metadata_samples_to_remove_as_string:
                    remove_file(f"{get_local_folder_for(settings.generated_dir_name, FileType.SequenceOrSampleData)}{x}.xml")
        except:
            logger.exception("Removal of metadata leftovers in the DB and XML files of the samples that failed was not successful.")

        database.try_py_function(vcm.insert_data_update_pipeline_event, pipeline_event)

# !!!!! ARE YOU LOOKING FOR prepared_parameters ? Use instead data_source.ncbi_any_virus.settings -> known_settings !!!!
#########################################################################################################################