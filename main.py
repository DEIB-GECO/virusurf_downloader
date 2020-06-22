import sys
from typing import Callable, Optional
import concurrent.futures
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
# logger.add(sink=sys.stderr,
#            level='TRACE',
#            format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
#                   "<level>{level: <8}</level> | "
#                   "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
#            colorize=True,
#            backtrace=True,
#            diagnose=True)
logger.add(sink=lambda msg: tqdm.write(msg, end=''),
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
database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=False)

#   ###################################     VIRUSES TO IMPORT    ###############
virus_taxon_ids = {
    'sars_cov2': 2697049
}
import_virus_sequences_in_parallel = True


def run():

    def import_virus_taxonomy_data(session: Session, virus_taxonomy_xml):
        virus_db_record = vcm.create_or_get_virus(session, virus_taxonomy_xml)
        return virus_db_record.virus_id

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
            nucleotide_variants = vcm.create_or_get_nucleotide_variants(session, sample, sequence, aligner)

            # TODO what about 'annotations' in ILcodiCE.parse_annotated_variants ?

    def try_import_virus_sequence(seq_acc_id):
        try:
            database_tom.try_py_function(import_virus_sequence, seq_acc_id)
        except:
            logger.exception(f'exception occurred while working on virus sample {seq_acc_id}.xml')

    for virus in virus_taxon_ids:
        # IMPORT VIRUS TAXON DATA
        logger.info(f'importing virus {virus}')
        virus_taxonomy_as_xml = vcm.download_virus_taxonomy_as_xml(virus_taxon_ids[virus])
        virus_id = database_tom.try_py_function(import_virus_taxonomy_data, virus_taxonomy_as_xml)

        # FIND VIRUS SEQUENCES
        logger.info(f'getting accession ids for virus sequences')
        refseq_accession_id, non_refseq_accession_ids = vcm.get_virus_sample_accession_ids(virus_taxon_ids[virus])
        # logger.warning('Two sequence accession ids are hardcoded in lines 63,64 of main.py. Uncomment line 62 to download '
        #                'all the sequence accession ids of this virus from NCBI.')
        # refseq_accession_id = 1798174254      # hardcoded value for tests
        # non_refseq_accession_ids = [1859094271]#, 1800242657, 1800242655, 1858732922, 1858732909]
        sequence_accession_ids = [refseq_accession_id] + non_refseq_accession_ids

        aligner: Optional[Callable] = None

        # # IMPORT VIRUS SEQUENCES
        logger.info(f'importing virus sequences and related tables')
        if import_virus_sequences_in_parallel:
            with tqdm(total=len(sequence_accession_ids)) as pbar:
                # import reference sequence first because it creates the nucleotide variant aligner
                try_import_virus_sequence(refseq_accession_id)
                pbar.update()
                # import other sequences
                with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                    future_generator = executor.map(try_import_virus_sequence, non_refseq_accession_ids)
                    for future_completed in future_generator:
                        pbar.update()
        else:
            for virus_seq_acc_id in tqdm(sequence_accession_ids):
                try_import_virus_sequence(virus_seq_acc_id)


run()
