import pickle
from collections import Counter, OrderedDict
from datetime import datetime
from decimal import Decimal
from typing import Tuple, Optional
from Bio import Entrez
from tqdm import tqdm
import data_cleaning_module
import stats_module
from pipeline_nuc_variants__annotations__aa import sequence_aligner
from loguru import logger
from time import sleep
from lxml import html, etree
import requests
from data_sources.ncbi_any_virus.ncbi_importer import AnyNCBIVNucSample
from xml_helper import text_at_node
from locations import get_local_folder_for, FileType
import wget
import os
from os.path import sep
from data_sources.ncbi_services import host_taxon_id_from_ncbi_taxon_name, download_ncbi_taxonomy_as_xml_from_name
from vcm import vcm
import dateutil.parser as dateparser
from geo_groups import geo_groups
from db_config import read_db_import_configuration as import_config, database_tom

Entrez.email = "example@mail.com"   # just to silence the warning. Then a correct email can be set later


cached_taxonomy = {}


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

    @staticmethod
    def is_positive_stranded():
        return True


class NMDCVirusSample:
    metadata_base_url = 'http://nmdc.cn/nCov/globalgenesequence/detail/'
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)

    def __init__(self, file_name: str):
        self._id = file_name.rstrip('.fasta')
        self.metadata = self._download_meta()
        self.nuc_seq = None     # cache for nucleotide sequence

    def __str__(self):
        return str(self.internal_id())

    def _download_meta(self):
        param = {
            "type": "globalgenesequence",
            "locus": self._id,
            "isnCov2": False,
            "id": ""
        }

        def make_request(attempts, error: Optional[Exception] = None):
            if attempts > 0:
                try:
                    res = requests.post('http://nmdc.cn/api/services/nmdcweb/api/search/detail', json=param)
                    if res.status_code != 200:
                        return make_request(attempts-1, ValueError(f"Fetch metadata response code is != 200. Full response: {res}"))
                    else:
                        return res.json()
                except Exception as e:
                    return make_request(attempts-1, e)
            else:
                logger.error(f'All attempts to fetch metadata for sample {self.internal_id()} failed.')
                raise error or Exception("")

        response = make_request(3)
        try:
            return response['data'][0]
        except IndexError:
            logger.error(f'Metadata for sample {self.internal_id()} not available')
            raise FileNotFoundError(f'Metadata for sample {self.internal_id()} not available')

    def internal_id(self):
        """
        An id used internally by this program when needed to log some errors. Choose whatever form and type you like.
        This id is not going to be written in the output data and possibly it will be printed only on the console in order
        to trace errors.
        """
        return self._id

    def primary_accession_number(self):
        return self._id

    def alternative_accession_number(self):
        return self.gisa_id()

    def strain(self):
        return self.metadata.get('isolate')

    @staticmethod
    def is_reference():
        return False

    def is_complete(self):
        def find_keyword(input_string: str):
            if input_string:
                lowered_string = input_string.lower()
                if 'complete' in lowered_string:
                    return True
                elif 'partial' in lowered_string:
                    return False
            else:
                return None

        res = find_keyword(self.metadata.get('definition'))
        if res is not None:
            return res
        res = find_keyword(self.metadata.get('description'))
        if res is not None:
            return res
        elif self.taxon_name() == 'Severe acute respiratory syndrome coronavirus 2' and self.length() < refseq_sc2_len*0.95:
            return False
        elif self.taxon_name() == 'Bat SARS-related coronavirus' and self.length() < refseq_sc1_len*0.95:
            return False
        else:
            return None

    def nucleotide_sequence(self):
        if not self.nuc_seq:
            try:
                with open(f'{fasta_folder}{self._id}.fasta', mode='r') as data:
                    data.readline()
                    self.nuc_seq = data.read().strip().replace('\n', '').lower()
            except Exception as e:
                logger.error(e)
                raise e
        return self.nuc_seq

    @staticmethod
    def strand():
        return 'positive'

    def length(self):
        length = self.metadata.get('glength') or len(self.nucleotide_sequence())
        return int(length)

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
        seq_tech = self.metadata.get('sequencingMethods')
        if seq_tech:
            return seq_tech.lower()
        else:
            return None

    def assembly_method(self):
        return self.metadata.get("jointMethods")

    @staticmethod
    def coverage():
        return None

    def collection_date(self):
        collection_date_format = self.metadata.get('collectionDateFormat')
        # logger.debug(f'collectionDateFormat: {collection_date_format}. collectionDate: {self.metadata.get("collectionDate")}. originalCollectionDate: {self.metadata.get("originalCollectionDate")}')
        if collection_date_format < '1900':     # for example, some unknown dates are encoded as '0001-01-01'
            return None
        elif collection_date_format:
            try:
                collection_date_format = dateparser.parse(collection_date_format, default=self.default_datetime).strftime('%Y-%m-%d')
            except dateparser._parser.ParserError as e:
                logger.warning(f'collection date string"{collection_date_format}" was parsed as {None}')
                return None
            return collection_date_format
        else:
            return None

    def isolation_source(self):
        return self.metadata.get('isolationSource')

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        input_string: str = self.metadata.get("country")
        geo_group = None
        country = None
        region = None
        if input_string:
            input_string.strip()
            if '/' in input_string:
                input_parts = input_string.split('/')
                input_parts = [i.strip() for i in input_parts]
                geo_group = input_parts[0]
                country = input_parts[1] if len(input_parts) > 1 else None
                region = input_parts[2] if len(input_parts) > 2 else None
            elif input_string.lower() in geo_groups.keys():
                country = input_string
                geo_group = geo_groups.get(country.lower())
            elif input_string == self.primary_accession_number():
                pass
            else:
                region = input_string
        return country, region, geo_group

    def submission_date(self):
        submit_date_format = self.metadata.get('submitDateFormat')
        # logger.debug(f'collectionDateFormat: {submit_date_format}. collectionDate: {self.metadata.get("submitDate")}.')
        if submit_date_format < '1900':     # for example, some unknown dates are encoded as '0001-01-01'
            return None
        elif submit_date_format:
            try:
                return datetime.strptime(submit_date_format, '%Y-%m-%d')
            except TypeError:
                return None
        else:
            return None

    def originating_lab(self):
        sampling_place: str = self.metadata.get('samplingPlace')
        if sampling_place:
            return sampling_place.replace(':', ' / ')
        else:
            return None

    def sequencing_lab(self) -> Optional[str]:
        return self.metadata.get('dept')

    def host_taxon_name(self) -> Optional[str]:
        name = self.metadata.get('host')
        if name is not None:
            name = name.strip()
            name = data_cleaning_module.correct_typos(name)
        return name

    def host_taxon_id(self) -> Optional[int]:
        try:
            return host_taxon_id_from_ncbi_taxon_name(self.host_taxon_name())
        except:
            logger.exception(f'Exception occurred while fetching the NCBI taxon id for taxon name {self.host_taxon_name()} in sample {self.internal_id()}')
            return None

    @staticmethod
    def gender() -> Optional[str]:
        return None

    @staticmethod
    def age() -> Optional[int]:
        return None

    @staticmethod
    def database_source() -> str:
        return 'NMDC'

    @staticmethod
    def bioproject_id():
        return None

    def taxon_name(self):
        name = self.metadata.get("spciesname")
        if name and name.lower() == 'human':
            name = 'Homo sapiens'
        return name

    def gisa_id(self):
        return self.metadata.get('gisaid')


def get_fasta_list():
    get_list_attempts = 3
    page = None
    while get_list_attempts > 0:
        try:
            page = requests.get(base_url)
            break
        except Exception as e:
            if get_list_attempts == 1:
                logger.exception('Attempt to fetch the list of fasta files continues to fail.')
                raise e
            logger.error('Attept to fetch the list of fasta files failed. New attmept in 10 seconds.')
            sleep(10)
            get_list_attempts -= 1
    tree = html.fromstring(page.content)
    file_list = tree.xpath('/html/body/pre/a')
    file_list = file_list[1:]     # removes '../'
    file_list = [text_at_node(a_node, '.', mandatory=True) for a_node in file_list]
    file_list = [item for item in file_list if item.startswith('NMDC')]       # keep only the sequences of SARS-COV
    return file_list


def download_fastas():
    def try_get_fasta_file(attempts: int):
        file_path = f'{fasta_folder}{file_name}'
        if not os.path.isfile(file_path):
            try:
                wget.download(f'http://nmdc.cn/SProject/virus/genbankfasta/{file_name}', out=file_path)
            except requests.exceptions.ChunkedEncodingError:
                if attempts > 0:
                    logger.error(f'Attempt to fetch of fasta file {file_name}. New attempt in 10 seconds.')
                    sleep(10)
                    try_get_fasta_file(attempts-1)
                else:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                    logger.exception(f'Download of fasta file {file_name} failed again. This file is being skipped.')
    for file_name in tqdm(fasta_list):
        try_get_fasta_file(3)


def reference_sequence(nuccore_query) -> str:
    # import a reference sequence from a different dataset that we'll use to call nucleotide variants
    # get reference sample accession id
    handle = Entrez.esearch(db="nuccore",
                            term=nuccore_query)
    response = Entrez.read(handle)
    handle.close()
    # Example response
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
    assert int(response['Count']) == 1, \
        "no reference sample found or multiple RefSeqs" + "please check: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch"
    reference_seq_id = response['IdList'][0]

    reference_seq_file_path = f"{get_local_folder_for(source_name='NMDC', _type=FileType.SequenceOrSampleData)}reference_sample.xml"
    if not os.path.exists(reference_seq_file_path):
        with Entrez.efetch(db="nuccore", id=reference_seq_id, rettype="gbc", retmode="xml") as handle:
            with open(reference_seq_file_path, 'w') as f:
                f.write(handle.read())

    reference_sample = AnyNCBIVNucSample(reference_seq_file_path, reference_seq_id)
    return reference_sample.nucleotide_sequence()


# GET FASTA LIST
base_url = 'http://nmdc.cn/SProject/virus/genbankfasta/'
fasta_list = []
# DOWNLOAD FASTA
fasta_folder: Optional[str] = None      # initialized elsewhere
taxonomy_folder: Optional[str] = None     # initialized elsewhere

# REFRENCE SEQUENCE
refseq_sc2 = None
refseq_sc2_len = 0
refseq_sc1 = None
refseq_sc1_len = 0

imported_viruses = set()


def import_samples_into_vcm():
    global fasta_list, refseq_sc1, refseq_sc2, refseq_sc1_len, refseq_sc2_len, cached_taxonomy, fasta_folder, \
        taxonomy_folder, imported_viruses
    db_params: dict = import_config.get_database_config_params()
    database_tom.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])
    fasta_folder = get_local_folder_for('NMDC', FileType.SequenceOrSampleData)
    taxonomy_folder = get_local_folder_for('NMDC', FileType.TaxonomyData)
    fasta_list = get_fasta_list()
    stats_module.schedule_samples(
        set_info=stats_module.StatsBasedOnIds([x.rstrip('.fasta') for x in fasta_list], True))
    logger.warning(
        f'{len(fasta_list)} files found at {base_url}. Some of them may be skipped because they have not metadata'
        f' or because they are not realted to SARS or SARS-coV2.')
    logger.info(f'downloading fasta files')
    download_fastas()

    # REFRENCE SEQUENCE
    refseq_sc2 = reference_sequence("(txid2697049[Organism]) AND srcdb_refseq[Properties]")
    refseq_sc2_len = len(refseq_sc2)
    refseq_sc1 = reference_sequence("txid694009[Organism:exp] NOT txid2697049[Organism] AND srcdb_refseq[Properties]")
    refseq_sc1_len = len(refseq_sc1)

    def virus_taxonomy_pipeline(session: database_tom, taxon: AnyNCBITaxon):
        return vcm.create_or_get_virus(session, taxon)

    # noinspection PyTypeChecker
    def metadata_pipeline(session: database_tom.Session, a_sample: NMDCVirusSample):
        try:
            experiment_id = vcm.create_or_get_experiment(session, a_sample)
            host_specie_id = vcm.create_or_get_host_specie(session, a_sample)
            host_sample_id = vcm.create_or_get_host_sample(session, a_sample, host_specie_id)
            sequencing_project_id = vcm.create_or_get_sequencing_project(session, a_sample)
            sequence = vcm.create_and_get_sequence(session, a_sample, virus_id, experiment_id, host_sample_id, sequencing_project_id)
            return sequence.sequence_id
        except Exception as e:
            if str(e).startswith('duplicate key value violates unique constraint "sequence_accession_id_key"'):
                logger.error(f'exception occurred while working on virus sample {a_sample}: {str(e)}')
            else:
                logger.exception(f'exception occurred while working on virus sample {a_sample}')
            raise database_tom.Rollback()

    def nucleotide__annotations__pipeline(session: database_tom.Session, a_sample: NMDCVirusSample, db_sequence_id):
        if a_sample.taxon_name() == 'Severe acute respiratory syndrome coronavirus 2':
            refseq = refseq_sc2
        elif a_sample.taxon_name() == 'Bat SARS-related coronavirus':
            refseq = refseq_sc1
        else:
            raise Exception(f'unknown taxon organism {a_sample.taxon_name()}')
        try:
            file_path = get_local_folder_for('NMDC', FileType.Annotations)+str(sample.primary_accession_number())+".pickle"
            if not os.path.exists(file_path):
                annotations_and_nuc_variants = sequence_aligner(
                    db_sequence_id,
                    refseq,
                    a_sample.nucleotide_sequence(),
                    'NC_045512',
                    f'.{sep}annotations{sep}new_ncbi_sars_cov_2.tsv',
                    'new_ncbi_sars_cov_2')
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
            stats_module.completed_sample(sample.primary_accession_number())
        except Exception:
            logger.exception(
                f'exception occurred while working on annotations and nuc_variants of virus sample '
                f'{a_sample.primary_accession_number()}. Rollback transaction.')
            raise database_tom.Rollback()

    logger.info('begin import of selected records')
    total_sequences_imported = 0
    total_sequences_skipped = 0
    log_of_gisaid_id_path = f"{get_local_folder_for('NMDC', FileType.Logs)}{os.path.sep}gisa_ids.txt"
    with open(log_of_gisaid_id_path, mode='w') as log_of_gisaid_id:
        for file in tqdm(fasta_list):
            try:
                sample = NMDCVirusSample(file)

                # filter samples by organism
                organism_name = sample.taxon_name()
                if organism_name != 'Severe acute respiratory syndrome coronavirus 2':
                    logger.info(f'Sample {file} skipped because related to organims {organism_name}')
                    total_sequences_skipped += 1
                    continue
                # download taxonomy for new organisms
                organism = cached_taxonomy.get('organism_name')
                if not organism:
                    organism_file = download_ncbi_taxonomy_as_xml_from_name(taxonomy_folder, organism_name)
                    organism = AnyNCBITaxon(organism_file)
                    cached_taxonomy[organism_name] = organism
            except FileNotFoundError:
                logger.error(f'Sample {file} skipped')
                total_sequences_skipped += 1
                continue
            except AssertionError:
                logger.exception(f'Sample {file} skipped')
                total_sequences_skipped += 1
                continue

            # virus id associated to this sample
            virus_id = database_tom.try_py_function(virus_taxonomy_pipeline, organism)
            if virus_id not in imported_viruses:
                imported_viruses.add(virus_id)
                database_tom.try_py_function(vcm.update_db_metadata, virus_id, 'NMDC')
            if virus_id:
                gisa_id = sample.gisa_id()
                if gisa_id:
                    log_of_gisaid_id.write(gisa_id+'\n')

                # import sample
                sequence_id = database_tom.try_py_function(metadata_pipeline, sample)
                if sequence_id:
                    database_tom.try_py_function(nucleotide__annotations__pipeline, sample, sequence_id)
                total_sequences_imported += 1

        logger.info(f'{total_sequences_imported} sequences imported.')
        logger.info(f'{total_sequences_skipped} sequences skipped.')
        logger.info(f'list of sequences with GISAID references at path: '+log_of_gisaid_id_path)

