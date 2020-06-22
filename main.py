import sys
import locations
from loguru import logger
import database_tom
import vcm as vcm
from virus_sample import VirusSample, download_or_get_virus_sample_as_xml
from sqlalchemy.orm.session import Session
from Bio import Entrez
import IlCodiCE
from tqdm import tqdm
Entrez.email = "Your.Name.Here@example.org"


#   #################################       PROGRAM ARGUMENTS   ##########################
wrong_arguments_message = 'The module main.py expects the following arguments:' \
                          'db_name, db_user', 'db_password', 'db_port.'
# noinspection PyBroadException
try:
    db_name = sys.argv[1]
    db_user = sys.argv[2]
    db_password = sys.argv[3]
    db_port = sys.argv[4]
except Exception:
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
logger.remove()  # removes default logger to stderr with level DEBUG
# choose what to print on console
logger.add(sys.stderr,
           level='TRACE',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True)

#   ###################################     SETUP FOLDERS   ###############################
locations.create_local_folders()

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
# init database
database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=True)


def run():

    # add a virus
    sars_cov2_taxonomy_id = 2697049  # FOR NOW USE JUST SARS-COV2 VIRUS
    logger.info(f'importing virus with taxonomy id {sars_cov2_taxonomy_id}')
    virus_taxonomy_as_xml = vcm.download_virus_taxonomy_as_xml(sars_cov2_taxonomy_id)

    def import_virus(session: Session):
        virus = vcm.create_or_get_virus(session, virus_taxonomy_as_xml)
        return virus.virus_id
    virus_id = database_tom.try_py_function(import_virus)

    # get associated sequence accession ids
    logger.info(f'getting accession ids for virus sequences')
    refseq_accession_id, non_refseq_accession_ids = vcm.get_virus_sample_accession_ids(sars_cov2_taxonomy_id)
    # logger.warning('Two sequence accession ids are hardcoded in lines 63,64 of main.py. Uncomment line 62 to download '
    #                'all the sequence accession ids of this virus from NCBI.')
    # refseq_accession_id = 1798174254      # hardcoded value for tests
    # non_refseq_accession_ids = [1829777589]
    sequence_accession_ids = [refseq_accession_id] + non_refseq_accession_ids
    aligner = None

    # import each sequence into the VCM
    def import_virus_sequence(session: Session, virus_seq_accession_id):
        virus_sample_file_path = download_or_get_virus_sample_as_xml(virus_seq_accession_id)
        sample = VirusSample(virus_sample_file_path, virus_seq_accession_id)

        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, virus_id, experiment, host_sample, sequencing_project)
        annotation = vcm.create_or_get_annotation(session, sample, sequence)

        nonlocal aligner
        if aligner is None and sample.is_reference():
            aligner = IlCodiCE.create_aligner_to_reference(reference=sequence.nucleotide_sequence, annotation_file='sars_cov_2_annotations.tsv', is_gisaid=False)
            logger.info('sequence aligner generated')
        else:
            # logger.info(f'doing global alignment of sequence accession id {virus_seq_accession_id}')

            # annotated_variants = aligner(sequence=sequence.nucleotide_sequence, sequence_id=virus_seq_accession_id)
            # print('annotated variants')
            # print(annotated_variants)

            # TODO row generated from parse_annotated_variants is incompatible with the table schema of Table 'variant'.
            # TODO fix table schema or function
            # nucleotide_variants = IlCodiCE.parse_annotated_variants(aligner)
            # print('nucleotide variants')
            # print(nucleotide_variants)
            # sys.exit(0)
            pass

    logger.info(f'importing virus sequences and related tables')
    for virus_seq_acc_id in tqdm(sequence_accession_ids):
        try:
            database_tom.try_py_function(import_virus_sequence, virus_seq_acc_id)
        except:
            logger.exception(f'exception occurred while working on virus sample {virus_seq_acc_id}.xml')


run()
