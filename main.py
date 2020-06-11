import sys
import locations
from loguru import logger
import vcm
from virus_sample import VirusSample, download_virus_sample_as_xml
from Bio import Entrez
Entrez.email = "Your.Name.Here@example.org"

#   ###################################      SETUP LOGGER    ##############################
logger.remove()  # removes default logger to stderr with level DEBUG
# choose what to print on console
logger.add(sys.stderr,
           level='TRACE',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<blue>{extra[request_id]}</blue> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True)

#   ###################################     SETUP FOLDERS   ###############################
locations.create_local_folders()

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
sars_cov2_taxonomy_id = 2697049   # FOR NOW USE JUST SARS-COV2 VIRUS
virus_taxonomy_as_xml = vcm.download_virus_taxonomy_as_xml(sars_cov2_taxonomy_id)
virus = vcm.create_virus(virus_taxonomy_as_xml)

refseq_accession_id, non_refseq_accession_ids = vcm.get_virus_sample_accession_ids(sars_cov2_taxonomy_id)
# refseq_accession_id = 1798174254      # hardcoded value for tests

virus_sample_file_path = download_virus_sample_as_xml(refseq_accession_id)
sample = VirusSample(virus_sample_file_path, refseq_accession_id)

experiment = vcm.create_or_get_experiment(sample)
host_sample = vcm.create_or_get_host_sample(sample)
sequencing_project = vcm.create_or_get_sequencing_project(sample)
sequence = vcm.create_sequence(sample, virus, experiment, host_sample, sequencing_project)
annotation = vcm.create_or_get_annotation(sample, sequence)
# TODO nucleotide_variant
# TODO amino acid variant

