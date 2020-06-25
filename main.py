import sys

from sqlalchemy.exc import SQLAlchemyError
from typing import Callable, Optional, List, Tuple
import locations
from loguru import logger
import database_tom
import vcm as vcm
from virus_sample import VirusSample, download_or_get_virus_sample_as_xml, DoNotImportSample, delete_virus_sample_xml
from multiprocessing import Queue, JoinableQueue, cpu_count, Process
import os
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
# on console print from level INFO on
logger.add(sink=lambda msg: tqdm.write(msg, end=''),
           level='INFO',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True)
# log to file any message of any security level
logger.add("./logs/log_{time}.log",
           level='TRACE',
           rotation='100 MB',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=False,
           backtrace=True,
           diagnose=True)

#   ###################################     SETUP FOLDERS   ###############################
locations.create_local_folders()

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
# init database
database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=True)

#   ###################################     VIRUSES TO IMPORT    ###############
virus_taxon_ids = {
    'sars_cov2': 2697049
}


class Sequential:

    def __init__(self, virus_id: int):
        self.virus_id = virus_id
        self.aligner: Optional[Callable] = None
        self.reference_sample: Optional[VirusSample] = None

    def import_virus_sample(self, session: Session, virus_seq_accession_id):
        virus_sample_file_path = download_or_get_virus_sample_as_xml(virus_seq_accession_id)
        sample = VirusSample(virus_sample_file_path, virus_seq_accession_id)

        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, self.virus_id, experiment, host_sample,
                                              sequencing_project)
        vcm.create_annotation_and_aa_variants(session, sample, sequence, self.reference_sample)

        if self.aligner is None and sample.is_reference():
            self.aligner = IlCodiCE.create_aligner_to_reference(reference=sequence.nucleotide_sequence,
                                                                annotation_file='sars_cov_2_annotations.tsv',
                                                                is_gisaid=False)
            self.reference_sample = sample
            logger.info('reference sequence imported')
        else:
            vcm.create_nucleotide_variants_and_impacts(session, sample, sequence.sequence_id, self.aligner)


class Parallel:

    MAX_PROCESSES = 5
    MAX_QUEUE = 6

    def __init__(self, virus_id: int):
        self.virus_id = virus_id
        self.reference_sample: Optional[VirusSample] = None
        # empty job queue
        self._queue = JoinableQueue(Parallel.MAX_QUEUE)
        self.workers: List[Parallel.Consumer] = []
        self.shared_session: Optional[Session] = None

    def import_virus_sample(self, session: Session, virus_seq_accession_id):
        # do this synchronously
        virus_sample_file_path = download_or_get_virus_sample_as_xml(virus_seq_accession_id)
        sample = VirusSample(virus_sample_file_path, virus_seq_accession_id)

        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, self.virus_id, experiment, host_sample,
                                              sequencing_project)
        vcm.create_annotation_and_aa_variants(session, sample, sequence, self.reference_sample)

        if not self.workers and sample.is_reference():
            self.reference_sample = sample
            self.shared_session = database_tom.get_session()
            # prepare workers
            for i in range(Parallel.number_of_processes()):
                aligner = IlCodiCE.create_aligner_to_reference(reference=sequence.nucleotide_sequence,
                                                               annotation_file='sars_cov_2_annotations.tsv',
                                                               is_gisaid=False)
                worker = Parallel.Consumer(self._queue, aligner, self.shared_session)
                self.workers.append(worker)
                worker.start()
            logger.info('reference sequence imported')
        else:
            # do this asynchronously
            sample.cache_nucleotide_sequence_and_release_source_xml_reference()  # release xml reference to prevent errors with multiprocessing
            logger.info(f'queue size {self._queue.qsize()}')
            self._queue.put([sample, sequence.sequence_id])
            logger.info(f'nuc variants for sample {virus_seq_accession_id} scheduled')

    class Consumer(Process):
        def __init__(self, jobs: JoinableQueue, aligner: Callable, shared_session: Session):
            super().__init__()
            self.jobs = jobs
            self.aligner: Optional[Callable] = aligner
            self.shared_session: Session = shared_session
            logger.info('worker started')

        def run(self):
            while True:
                element = self.jobs.get(block=True, timeout=None)
                if element is None:
                    # shutdown this worker
                    print('shutting down a consumer')
                    self.jobs.task_done()
                    break
                else:
                    try:
                        vcm.create_nucleotide_variants_and_impacts(self.shared_session, element[0], element[1], self.aligner)
                        self.shared_session.commit()
                        logger.info(f'nucleotide variants for sample {element[0].alternative_accession_number()}')
                    except SQLAlchemyError as e:
                        self.shared_session.rollback()
                    self.jobs.task_done()

    @staticmethod
    def number_of_processes():
        max_p = []
        # get max allowed processes
        try:
            max_p.append(len(os.sched_getaffinity(0)))
        except AttributeError:
            pass    # method not available on this system (can happen on Windows and macOS)
        # get number of CPUs
        try:
            max_p.append(cpu_count())
        except NotImplementedError:
            pass
        # take the smallest number of above two and Parallel.MAX_PROCESSES
        max_p.append(Parallel.MAX_PROCESSES)
        used_processes = min(max_p)
        if used_processes < Parallel.MAX_PROCESSES:
            logger.warning(f'maximum number of CPUs or usable process reached. At most {used_processes} will be used.')
        return used_processes

    def tear_down(self):
        # terminate workers
        for w in self.workers:
            self._queue.put(None, block=True, timeout=None)
        logger.info('waiting all workers to finish')
        self._queue.join()
        logger.info('all processes_finished')
        self.shared_session.close()


def run():

    def import_virus_taxonomy_data(session: Session, virus_taxonomy_xml):
        virus_db_record = vcm.create_or_get_virus(session, virus_taxonomy_xml)
        return virus_db_record.virus_id

    def try_import_virus_sample(seq_acc_id):
        nonlocal successful_imports
        try:
            database_tom.try_py_function(import_method.import_virus_sample, seq_acc_id)
            successful_imports += 1
        except DoNotImportSample:
            delete_virus_sample_xml(seq_acc_id)
            logger.info(f'virus sample with accession id {seq_acc_id} has not been imported and was removed from disk')
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
        # logger.warning('Sequence accession ids are hardcoded in main.py.')
        # refseq_accession_id = 1798174254      # hardcoded value for tests
        # non_refseq_accession_ids = [1852393386]#, 1852393360, 1852393373,  1852393398, 1859035944, 1859035892, 1800242649, 1800242657, 1858732896, 1799706760, 1800242651, 1800242659, 1858732909, 1799706762, 1800242653, 1800242661, 1858732922, 1800242639, 1800242655, 1850952215, 1859094271]
        non_refseq_accession_ids = non_refseq_accession_ids[::-1]
        sequence_accession_ids = [refseq_accession_id] + non_refseq_accession_ids

        # # IMPORT VIRUS SEQUENCES
        logger.info(f'importing virus sequences and related tables')
        global import_method
        import_method = import_method(virus_id)

        successful_imports = 0
        for virus_seq_acc_id in tqdm(sequence_accession_ids):
            try_import_virus_sample(virus_seq_acc_id)

        logger.info('main process completed')
        import_method.tear_down()
        logger.info(f'successful imports: {successful_imports} (not reliable when parallel processing)')


import_method = Parallel
run()
