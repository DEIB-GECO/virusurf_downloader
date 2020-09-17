import os
from collections import OrderedDict
from data_sources.ncbi_any_virus.ncbi_importer import AnyNCBIVNucSample

from lxml import etree
import wget as wget
from Bio import Entrez
from tqdm import tqdm
from typing import Generator, Callable, Optional
from loguru import logger
from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from locations import get_local_folder_for, FileType
import database_tom
from vcm import sequence_primary_accession_ids, remove_sequence_and_meta
import stats_module
from xml_helper import text_at_node


class COGUKSarsCov2:

    name = 'COG-UK_sars_cov_2'
    sequence_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-09-03/cog_2020-09-03_sequences.fasta'
    metadata_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-09-03/cog_2020-09-03_metadata.csv'

    def __init__(self):
        logger.info(f'importing virus {self.name} using NCBI SC2 for taxonomy data')
        self.tax_tree = download_virus_taxonomy_as_xml(
            get_local_folder_for(source_name=self.name, _type=FileType.TaxonomyData),
            self.taxon_id())
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
        id_local_samples = set(database_tom.try_py_function(
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
                        yield COGUKSarsCov2Sample(sample)

    def reference_sequence(self) -> str:
        if not hasattr(self, 'reference_sample'):
            # import a reference sequence from a different dataset that we'll use to call nucleotide variants
            # get reference sample accession id
            handle = Entrez.esearch(db="nuccore",
                                    term=f"(txid{self.taxon_id()}[Organism]) AND srcdb_refseq[Properties]")
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

            reference_seq_file_path = f"{get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)}reference_sample.xml"
            if not os.path.exists(reference_seq_file_path):
                with Entrez.efetch(db="nuccore", id=reference_seq_id, rettype="gbc", retmode="xml") as handle, open(
                        reference_seq_file_path, 'w') as f:
                    for line in handle:
                        f.write(line)

            self.reference_sample = AnyNCBIVNucSample(reference_seq_file_path, reference_seq_id)
        return self.reference_sample.nucleotide_sequence()


def download_or_get_sample_data(containing_directory: str) -> (str, str):
    """
    :return: the local file path of the downloaded sequence and metadata files.
    """
    sequence_local_file_path = containing_directory + COGUKSarsCov2.sequence_file_url.rsplit('/', maxsplit=1)[1]
    metadata_local_file_path = containing_directory + COGUKSarsCov2.metadata_file_url.rsplit('/', maxsplit=1)[1]
    if not os.path.exists(sequence_local_file_path) or not os.path.exists(metadata_local_file_path):
        logger.info(f'downloading sample sequences for COG-UK data from {COGUKSarsCov2.sequence_file_url} ...')
        wget.download(COGUKSarsCov2.sequence_file_url, sequence_local_file_path)
        logger.info(f'downloading sample metadata for COG-UK data from {COGUKSarsCov2.metadata_file_url} ...')
        wget.download(COGUKSarsCov2.metadata_file_url, metadata_local_file_path)
    return sequence_local_file_path, metadata_local_file_path


def delete_virus_sample_xml(containing_directory: str, sample_accession_id: int):
    """
    :param containing_directory: directory where the file resides
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    """
    local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
    try:
        os.remove(local_file_path)
    except OSError as e:
        logger.error(f"Failed to remove file {local_file_path} with error: {e.strerror}")


def download_or_get_virus_sample_as_xml(containing_directory: str, sample_accession_id: int) -> str:
    """
    :param containing_directory: directory where the file will be downloaded and cached
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
    with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml") as handle, open(local_file_path, 'w') as f:
        for line in handle:
            f.write(line)
    return local_file_path


def download_virus_taxonomy_as_xml(containing_directory: str, taxon_id) -> etree.ElementTree:
    # write taxonomy tree
    local_file_path = f"{containing_directory}{os.path.sep}{taxon_id}.xml"
    if not os.path.exists(local_file_path):
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                open(local_file_path, 'w') as f:
            for line in handle:
                f.write(line)
    tax_tree = etree.parse(local_file_path, parser=etree.XMLParser(remove_blank_text=True))
    return tax_tree
