import os
import sys

import wget as wget
from Bio import Entrez
from tqdm import tqdm
from typing import Generator, Callable

from loguru import logger

from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from data_sources.ncbi_sars_cov_2.sample import NCBISarsCov2Sample
from data_sources.ncbi_sars_cov_2.virus import NCBISarsCov2
from data_sources.ncbi_sars_cov_2.virus import download_or_get_virus_sample_as_xml, get_virus_sample_accession_ids
from locations import local_folder


class COGUKSarsCov2(NCBISarsCov2):

    name = 'COGUK_sars_cov_2'
    sequence_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_sequences.fasta'
    metadata_file_url = 'https://cog-uk.s3.climb.ac.uk/2020-05-08/cog_2020-05-08_metadata.csv'

    def __init__(self):
        logger.info(f'importing virus {COGUKSarsCov2.name} using {super().name} for taxonomy data')
        super().__init__()
        self.sequence_file_path, self.metadata_file_path = download_or_get_sample_data()

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

    def virus_samples(self) -> Generator[COGUKSarsCov2Sample, None, None]:
        # import metadata (few KB, we can save it in memory)
        # noinspection PyAttributeOutsideInit
        logger.info('parsing metadata file...')
        meta = {}
        with open(self.metadata_file_path, mode='r') as metadata_file:
            metadata_file.readline()    # skip header
            for line in metadata_file:
                try:
                    key, content = line.split(sep=',', maxsplit=1)
                    meta[key] = content.rstrip()
                    # logger.debug(f'new meta_key: {key}')
                except:
                    logger.error(f'Unable to parse the following line from the metadata file {self.metadata_file_path}:\n\t'
                                 f'{line}')

        # generate a VirusSample for each line from region data file
        logger.info('reading sample sequence file...')
        with open(self.sequence_file_path, mode='r') as seq_file:
            progress = tqdm(total=len(meta.keys()))
            while True:
                progress.update()
                sample_key = seq_file.readline()
                sample_sequence = seq_file.readline()
                if not sample_sequence or not sample_key:
                    break  # EOF
                else:
                    sample_key = sample_key[1:].rstrip()  # in the sequence file, the strain name is preceded by a '>', while in the metadata file not
                    sample = {
                        COGUKSarsCov2Sample.STRAIN_NAME: sample_key,
                        COGUKSarsCov2Sample.NUC_SEQUENCE: sample_sequence.rstrip()
                    }
                    try:
                        sample[COGUKSarsCov2Sample.METADATA_RAW_STRING] = meta[sample_key]
                    except KeyError:
                        logger.error(f'Found a sequence without a paired metadata string. Foreign key was {sample_key}')
                    yield COGUKSarsCov2Sample(sample)

    def nucleotide_variant_aligner(self) -> Callable:
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

            reference_seq_file_path = f"{local_folder}{os.path.sep}coguk{os.path.sep}reference_sample.xml"
            if not os.path.exists(reference_seq_file_path):
                with Entrez.efetch(db="nuccore", id=reference_seq_id, rettype="gbc", retmode="xml") as handle, open(
                        reference_seq_file_path, 'w') as f:
                    for line in handle:
                        f.write(line)

            self.reference_sample = NCBISarsCov2Sample(reference_seq_file_path, reference_seq_id)
        return self.reference_sample.nucleotide_var_aligner()


def download_or_get_sample_data() -> (str, str):
    """
    :return: the local file path of the downloaded sequence and metadata files.
    """
    directory = f"{local_folder}{os.sep}coguk{os.sep}"
    sequence_local_file_path = directory + COGUKSarsCov2.sequence_file_url.rsplit('/', maxsplit=1)[1]
    metadata_local_file_path = directory + COGUKSarsCov2.metadata_file_url.rsplit('/', maxsplit=1)[1]
    if not os.path.exists(sequence_local_file_path) or not os.path.exists(metadata_local_file_path):
        if not os.path.exists(directory):
            os.makedirs(directory)
        logger.info(f'downloading sample sequences for COG-UK data from {COGUKSarsCov2.sequence_file_url} ...')
        wget.download(COGUKSarsCov2.sequence_file_url, sequence_local_file_path)
        logger.info(f'downloading sample metadata for COG-UK data from {COGUKSarsCov2.metadata_file_url} ...')
        wget.download(COGUKSarsCov2.metadata_file_url, metadata_local_file_path)
    return sequence_local_file_path, metadata_local_file_path
