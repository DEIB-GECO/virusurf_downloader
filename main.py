import sys
import os
from loguru import logger
from typing import Callable, Optional, List
import database_tom
import vcm as vcm
from data_sources.ncbi_sars_cov_1.virus import NCBISarsCov1
from data_sources.virus_sample import VirusSample
from data_sources.virus import VirusSource
from data_sources.gisaid_sars_cov_2.virus import GISAIDSarsCov2
from data_sources.ncbi_sars_cov_2.virus import NCBISarsCov2
from data_sources.ncbi_any_virus.ncbi_importer import prepared_parameters, import_samples_into_vcm_except_annotations_nuc_vars
import data_sources.nmdc.procedure as nmdc
from multiprocessing import JoinableQueue, cpu_count, Process
from sqlalchemy.orm.session import Session
from Bio import Entrez
from tqdm import tqdm
import warnings

Entrez.email = "Your.Name.Here@example.org"


#   #################################       PROGRAM ARGUMENTS   ##########################
wrong_arguments_message = 'The module main.py expects the following arguments:' \
                          'db_name, recreate_db?, db_user, db_password, db_port, source_to_import\n' \
                          'source_to_import accept values:\n' \
                          'genbank-sars-cov2\n' \
                          'cog-uk\n' \
                          'gisaid\n' \
                          'genbank-sars-cov1\n'
for p in prepared_parameters.keys():
    wrong_arguments_message += f'{p}\n'
wrong_arguments_message += 'just_make_indexes\n'
wrong_arguments_message += 'optionally you can specify a range of samples to import as <min> <max>'
# noinspection PyBroadException
try:
    db_name = sys.argv[1]
    _db_recreate = sys.argv[2].lower()
    db_recreate = True if _db_recreate == 'yes' or _db_recreate == 'true' else False
    db_user = sys.argv[3]
    db_password = sys.argv[4]
    db_port = sys.argv[5]
    source = sys.argv[6].lower()
    try:
        _from = int(sys.argv[7])
        to = int(sys.argv[8])
    except Exception:
        _from = None
        to = None
except Exception:
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
logger.remove()  # removes default logger to stderr with level DEBUG
# on console print from level INFO on
logger.add(sink=lambda msg: tqdm.write(msg, end=''),
           level='TRACE',   #TODO change before committing
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True,
           enqueue=True)
# log to file any message of any security level
logger.add("./logs/log_"+source+"_{time}.log",
           level='TRACE',
           rotation='100 MB',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=False,
           backtrace=True,
           diagnose=True,
           enqueue=True)
# redirect warnings
def customwarn(message, category, filename, lineno, file=None, line=None):
    logger.warning(warnings.formatwarning(message, category, filename, lineno))
warnings.showwarning = customwarn

logger.info(f"main.py {' '.join(sys.argv[1:])}")

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
# init database
database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=db_recreate)

#   ###################################     VIRUSES TO IMPORT    ###############
viruses: List[VirusSource] = [

]


class Sequential:

    def __init__(self, virus_id: int):
        self.virus_id = virus_id
        self.aligner: Optional[Callable] = None
        self.reference_sample: Optional[VirusSample] = None

    def import_virus_sample(self, session: Session, sample: VirusSample):
        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, self.virus_id, experiment, host_sample,
                                              sequencing_project)
        vcm.create_annotation_and_aa_variants(session, sample, sequence, self.reference_sample)

        if self.aligner is None and sample.is_reference():
            self.aligner = sample.nucleotide_var_aligner()
            self.reference_sample = sample
            logger.info('reference sequence imported')
        else:
            vcm.create_nucleotide_variants_and_impacts(session, sample, sequence.sequence_id, self.aligner)

    def tear_down(self):
        pass

class Parallel:

    MAX_PROCESSES = 3

    def __init__(self, virus_id: int):
        self.virus_id = virus_id
        self.reference_sample: Optional[VirusSample] = None
        # empty job queue
        queue_size = self.number_of_processes()+10
        self._queue = JoinableQueue(maxsize=queue_size)
        logger.info(f'Queue size set to accept at most {queue_size} before pausing the producer.')
        self.workers: List[Parallel.Consumer] = []
        self.shared_session: Optional[Session] = None

    def import_virus_sample(self, session: Session, sample):
        # do this synchronously
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
            for _ in range(Parallel.number_of_processes()):
                worker = Parallel.Consumer(self._queue, sample.nucleotide_var_aligner(), self.shared_session)
                self.workers.append(worker)
                worker.start()
            logger.info('reference sequence imported')
        else:
            # schedule nucleotide variants to be called asynchronously
            sample.on_before_multiprocessing()
            self._queue.put([sample, sequence.sequence_id])
            logger.info(f'nucleotide variant calling for sequence {sample.internal_id()} scheduled')

    class Consumer(Process):
        def __init__(self, jobs: JoinableQueue, aligner: Callable, shared_session: Session):
            super().__init__()
            self.jobs = jobs
            self.aligner: Optional[Callable] = aligner
            self.shared_session: Session = shared_session
            logger.info('worker started')

        def run(self):
            while True:
                job = self.jobs.get(block=True, timeout=None)
                if job is None:
                    logger.info('shutting down a consumer')
                    self.jobs.task_done()
                    break
                else:
                    sample, sequence_id = job
                    try:
                        vcm.create_nucleotide_variants_and_impacts(self.shared_session, sample, sequence_id, self.aligner)
                        self.shared_session.commit()
                        logger.info(f'nucleotides of sequence {sample.internal_id()} imported')
                    except:
                        logger.exception(f'unknown exception while calling nucleotides on sample {sample.internal_id()}')
                        database_tom.rollback(self.shared_session)
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
        for _ in self.workers:
            self._queue.put(None, block=True, timeout=None)
        logger.info('waiting all workers to finish')
        self._queue.join()
        logger.info('all processes_finished')
        self.shared_session.close()


def run():

    def import_virus(session: Session, virus: VirusSource):
        return vcm.create_or_get_virus(session, virus).virus_id

    def try_import_virus_sample(sample: VirusSample):
        nonlocal successful_imports
        try:
            database_tom.try_py_function(import_method.import_virus_sample, sample)
            successful_imports += 1
        except:
            logger.exception(f'exception occurred while working on virus sample {sample.internal_id()}')

    for virus in viruses:
        # IMPORT VIRUS TAXON DATA
        virus_id = database_tom.try_py_function(import_virus, virus)

        # # IMPORT VIRUS SEQUENCES
        logger.info(f'importing virus sequences and related tables')
        global import_method
        import_method = import_method(virus_id)
        successful_imports = 0

        for s in virus.virus_samples():
            try_import_virus_sample(s)

        logger.info('main process completed')
        import_method.tear_down()
        logger.info(f'successful imports: {successful_imports} (not reliable when parallel processing)')

try:
    if 'index' in source:
        database_tom.create_indexes()
    elif source in ['coguk', 'cog-uk']:
        import data_sources.coguk_sars_cov_2.procedure
    elif source in ['ncbi-sars-cov2', 'genbank-sars-cov2', 'genbank-sarscov2', 'genbank-sars-cov-2']:
        import_method = Parallel
        viruses = [NCBISarsCov2()]
        run()
    elif source == 'gisaid':
        import_method = Sequential
        viruses = [GISAIDSarsCov2()]
        run()
    elif source in ['ncbi-sars-cov1', 'genbank-sars-cov1', 'genbank-sarscov1', 'genbank-sars-cov']:
        import_method = Parallel
        viruses = [NCBISarsCov1()]
        run()
    elif source in prepared_parameters.keys():
        import_samples_into_vcm_except_annotations_nuc_vars(*prepared_parameters[source], from_sample=_from, to_sample=to)
    elif 'nmdc' in source:
        nmdc.import_samples_into_vcm()
    else:
        logger.error(f'the argument {source} is not recognised.\n'+wrong_arguments_message)
except:
    logger.exception('FATAL ERROR')  # this is just to make sure the exception is written to the log file before crashing
    sys.exit(-1)
finally:
    logger.complete()
