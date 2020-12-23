import os
import shutil
from collections import OrderedDict
from data_sources.ncbi_any_virus.ncbi_importer import AnyNCBIVNucSample
from data_sources.ncbi_services import download_ncbi_taxonomy_as_xml, get_samples_accession_ids, download_or_get_ncbi_sample_as_xml
from data_sources.ncbi_any_virus.settings import known_settings as ncbi_known_settings
from lxml import etree
import wget as wget
from tqdm import tqdm
from typing import Generator, Optional
from loguru import logger
from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from locations import get_local_folder_for, FileType, remove_file
from db_config import database
from vcm.vcm import sequence_primary_accession_ids
import stats_module
from xml_helper import text_at_node
from urllib.request import Request, urlopen


class COGUKSarsCov2:

    name = 'coguk_sars_cov_2'
    sequence_file_url = 'https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_all.fasta'
    metadata_file_url = 'https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv'

    def __init__(self):
        logger.info(f'importing virus {self.name} using NCBI SC2 for taxonomy data')
        # fetch taxonomy data from NCBI
        taxonomy_file_path = download_ncbi_taxonomy_as_xml(
            get_local_folder_for(source_name=self.name, _type=FileType.TaxonomyData),
            self.taxon_id())
        try:
            self.tax_tree = etree.parse(taxonomy_file_path, parser=etree.XMLParser(remove_blank_text=True))
        except etree.XMLSyntaxError as e:  # happens on AWS if for some reason the downloaded file is corrupted
            remove_file(taxonomy_file_path)
            ncbi_sc2_taxonomy_dir = get_local_folder_for(
                source_name=ncbi_known_settings["sars_cov_2"]["generated_dir_name"], _type=FileType.TaxonomyData)
            alternative_taxonomy_path = ncbi_sc2_taxonomy_dir + f"{ncbi_known_settings['sars_cov_2']['virus_taxon_id']}.xml"
            if os.path.exists(alternative_taxonomy_path):
                shutil.copyfile(alternative_taxonomy_path, taxonomy_file_path)
                self.tax_tree = etree.parse(taxonomy_file_path, parser=etree.XMLParser(remove_blank_text=True))
            else:
                logger.error(f"Taxonomy file of SARS-CoV-2 was empty. Attempt to use the one from {ncbi_sc2_taxonomy_dir} "
                             f"failed because the filed doesn't exist. Can't proceed.")
                raise e
        # fetch latest source data
        download_dir = get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)
        self.sequence_file_path, self.metadata_file_path = download_or_get_sample_data(download_dir)

    def taxon_id(self):
        return 2697049

    def taxon_name(self):
        return text_at_node(self.tax_tree, './Taxon/ScientificName')

    def family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName')

    def sub_family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName')

    def genus(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName')

    def species(self):
        # species_taxon_id = text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/TaxId')
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName',
                            mandatory=False)

    def equivalent_names(self):
        genbank_acronym = text_at_node(self.tax_tree, './/GenbankAcronym', mandatory=False)
        equivalent_names = self.tax_tree.xpath('.//EquivalentName')
        equivalent_names = [x.text for x in equivalent_names]
        if genbank_acronym:
            equivalent_names.insert(0, genbank_acronym)
        equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
        equivalent_names = ", ".join(equivalent_names)
        return equivalent_names

    def molecule_type(self):
        return 'RNA'

    def is_single_stranded(self):
        return True

    def is_positive_stranded(self):
        return True

    @staticmethod
    def _deltas(virus_id, id_remote_samples):
        """
        :param virus_id: virus_id from the DB of the virus being imported
        :param id_remote_samples: alternative accession id of the sequences avbailable from remote
        :return: three sets containing the alternative accession ids of
            1. sequences that are unchanged (local == remote)
            2. outdated sequences (to be removed)
            3. sequences new and to be imported
        """
        id_local_samples = set(database.try_py_function(
            sequence_primary_accession_ids, virus_id, 'COG-UK'
        ))
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

    def virus_samples(self, virus_id: int, from_sample: Optional[int] = None, to_sample: Optional[int] = None) -> Generator[COGUKSarsCov2Sample, None, None]:
        # import metadata (few KB, we can save it in memory)
        # noinspection PyAttributeOutsideInit
        logger.info('parsing metadata file...')
        meta = {}
        with open(self.metadata_file_path, mode='r') as metadata_file:
            metadata_file.readline()    # skip header
            for line in metadata_file:
                try:
                    key, content = line.split(sep=',', maxsplit=1)  # key = <sequence name>, content = <other metadata>
                    meta[key] = content.rstrip()
                    # logger.debug(f'new meta_key: {key}')
                except:
                    logger.error(f'Unable to parse the following line from the metadata file {self.metadata_file_path}:\n\t'
                                 f'{line}')

        # do _deltas() or import everything
        id_new_sequences = meta.keys()

        # slice the sequences
        id_new_sequences = sorted(list(id_new_sequences))
        if from_sample is not None and to_sample is not None:
            id_new_sequences = id_new_sequences[from_sample:to_sample]

        # proceed importing selected sequences
        meta = dict(filter(lambda kv: kv[0] in id_new_sequences, meta.items()))

        stats_module.schedule_samples(
            stats_module.StatsBasedOnIds(id_new_sequences, True, virus_id, ['COG-UK']))

        # generate a VirusSample for each line from region data file
        logger.info('reading sample sequence file...')
        with open(self.sequence_file_path, mode='r') as seq_file:
            progress = tqdm(total=len(meta.keys()))
            while True:
                sample_key = seq_file.readline()
                sample_sequence = seq_file.readline()
                if not sample_sequence or not sample_key:
                    break  # EOF
                else:
                    sample_key = sample_key[1:].rstrip()  # in the sequence file, the strain name is preceded by a '>', while in the metadata file not
                    if sample_key not in meta.keys():
                        continue
                    else:
                        progress.update()
                        sample = {
                            COGUKSarsCov2Sample.STRAIN_NAME: sample_key,
                            COGUKSarsCov2Sample.NUC_SEQUENCE: sample_sequence.rstrip()
                        }
                        try:
                            sample[COGUKSarsCov2Sample.METADATA_RAW_STRING] = meta[sample_key]
                        except KeyError:
                            logger.error(
                                f'Found a sequence without a paired metadata string. Foreign key was {sample_key}')
                        try:
                            yield COGUKSarsCov2Sample(sample)
                        except:
                            logger.exception(f"Sample {sample_key} skipped due to an error while parsing the input data.")

    def reference_sequence(self) -> str:
        if not hasattr(self, 'reference_sample'):
            # import a reference sequence from a different dataset that we'll use to call nucleotide variants
            # get reference sample accession id
            ncbi_reference_sample_query = ncbi_known_settings["sars_cov_2"]["reference_sample_query"]
            reference_accession_id = get_samples_accession_ids(ncbi_reference_sample_query)
            assert len(reference_accession_id) == 1, \
                "no reference sample found or multiple RefSeqs. Please correct the query used on NCBI nuccore"
            # download file as XML
            reference_seq_file_path = download_or_get_ncbi_sample_as_xml(
                get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData),
                reference_accession_id[0]
            )
            # parse file and cache the object
            self.reference_sample = AnyNCBIVNucSample(reference_seq_file_path, reference_accession_id[0])
        return self.reference_sample.nucleotide_sequence()


def download_or_get_sample_data(containing_directory: str) -> (str, str):
    """
    :return: the local file path of the downloaded sequence and metadata files.
    """
    def get_download_size(url) -> Optional[int]:
        """
        Returns the size of the downloadable resource in bytes supported by the protocol
        of the downloadable resource; None otherwise.
        """
        req = Request(url=url, method='HEAD')
        f = urlopen(req)
        if int(f.status) == 200:
            return int(f.headers['Content-Length'])
        else:
            return None

    def get_local_file_size(path) -> Optional[int]:
        """
        Returns the size of the local file in bytes if it exists; None otherwise.
        """
        if os.path.exists(path):
            return os.stat(path).st_size
        else:
            return None

    def download_coguk_data():
        logger.info(f'downloading sample sequences for COG-UK data from {COGUKSarsCov2.sequence_file_url} ...')
        wget.download(COGUKSarsCov2.sequence_file_url, sequence_local_file_path)
        logger.info(f'\ndownloading sample metadata for COG-UK data from {COGUKSarsCov2.metadata_file_url} ...')
        wget.download(COGUKSarsCov2.metadata_file_url, metadata_local_file_path)
        logger.info('\n')

    sequence_local_file_path = containing_directory + COGUKSarsCov2.sequence_file_url.rsplit('/', maxsplit=1)[1]
    metadata_local_file_path = containing_directory + COGUKSarsCov2.metadata_file_url.rsplit('/', maxsplit=1)[1]
    if not os.path.exists(sequence_local_file_path) or not os.path.exists(metadata_local_file_path):
        download_coguk_data()
    else:
        # compare size of local with remote ones
        if get_download_size(COGUKSarsCov2.sequence_file_url) != get_local_file_size(sequence_local_file_path) or \
                get_download_size(COGUKSarsCov2.metadata_file_url) != get_local_file_size(metadata_local_file_path):
            download_coguk_data()
    return sequence_local_file_path, metadata_local_file_path
