import os
import sys

import wget as wget
from Bio import Entrez
from tqdm import tqdm
from typing import Generator, Callable, Optional

from loguru import logger

from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from data_sources.ncbi_sars_cov_2.sample import NCBISarsCov2Sample
from data_sources.ncbi_sars_cov_2.virus import NCBISarsCov2
from locations import get_local_folder_for, FileType
import database_tom
from vcm import sequence_primary_accession_ids, remove_sequence_and_meta


class COGUKSarsCov2(NCBISarsCov2):

    name = 'COG-UK_sars_cov_2'
    sequence_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-07-28/cog_2020-07-28_sequences.fasta'
    metadata_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-07-28/cog_2020-07-28_metadata.csv'

    def __init__(self):
        logger.info(f'importing virus {self.name} using {super().name} for taxonomy data')
        super().__init__()
        download_dir = get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)
        self.sequence_file_path, self.metadata_file_path = download_or_get_sample_data(download_dir)

    def taxon_id(self):
        return super().taxon_id()

    def taxon_name(self):
        return super().taxon_name()

    def family(self):
        return super().family()

    def sub_family(self):
        return super().sub_family()

    def genus(self):
        return super().genus()

    def species(self):
        return super().species()

    def equivalent_names(self):
        return super().equivalent_names()

    def molecule_type(self):
        return super().molecule_type()

    def is_single_stranded(self):
        return super().is_single_stranded()

    def is_positive_stranded(self):
        return super().is_positive_stranded()

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

        # compare data with already imported data
        id_current_sequences = meta.keys()
        id_previously_imported_sequences = set(database_tom.try_py_function(
            sequence_primary_accession_ids, virus_id, 'COG-UK'
        ))
        id_outdated_sequences = id_previously_imported_sequences - id_current_sequences
        id_new_sequences = id_current_sequences - id_previously_imported_sequences
        logger.info(f'\n# total current sequences from source: {len(id_current_sequences)}. Of which\n'
                    f'# {len(id_previously_imported_sequences)} imported in previous runs\n'
                    f'# {len(id_new_sequences)} new sequences\n'
                    f'# {len(id_outdated_sequences)} are outdated and must be removed from previous import')

        # remove outdated sequences
        for prim_acc_id in id_outdated_sequences:
            database_tom.try_py_function(
                remove_sequence_and_meta, prim_acc_id, None
            )

        # proceed importing only new sequences
        meta = dict(filter(lambda kv: kv[0] in id_new_sequences, meta.items()))

        # generate a VirusSample for each line from region data file
        logger.info('reading sample sequence file...')
        counter = 0
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
                        if from_sample is not None and to_sample is not None:
                            if counter < from_sample:
                                progress.update()
                                counter += 1
                                continue
                            elif counter >= to_sample:
                                break
                            else:
                                counter += 1
                                progress.update()
                                sample = {
                                    COGUKSarsCov2Sample.STRAIN_NAME: sample_key,
                                    COGUKSarsCov2Sample.NUC_SEQUENCE: sample_sequence.rstrip()
                                }
                                try:
                                    sample[COGUKSarsCov2Sample.METADATA_RAW_STRING] = meta[sample_key]
                                except KeyError:
                                    logger.error(f'Found a sequence without a paired metadata string. Foreign key was {sample_key}')
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

            self.reference_sample = NCBISarsCov2Sample(reference_seq_file_path, reference_seq_id)
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
