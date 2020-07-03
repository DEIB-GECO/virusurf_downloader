import sys

import os
from loguru import logger
from typing import Callable, Optional, List

import locations
import database_tom
import vcm as vcm
from data_sources.virus_sample import VirusSample
from data_sources.virus import VirusSource
from data_sources.gisaid_sars_cov_2.virus import GISAIDSarsCov2
from data_sources.ncbi_sars_cov_2.virus import NCBISarsCov2
from data_sources.coguk_sars_cov_2.virus import COGUKSarsCov2
from multiprocessing import JoinableQueue, cpu_count, Process
from sqlalchemy.orm.session import Session
from Bio import Entrez
from tqdm import tqdm
Entrez.email = "Your.Name.Here@example.org"


class Sequential:

    def import_virus_sample(self, session: Session, sample: VirusSample):
        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, virus_id, experiment, host_sample,
                                              sequencing_project)

        # vcm.create_nucleotide_variants_and_impacts(session, sample, sequence.sequence_id, virus.nucleotide_variant_aligner())

    def tear_down(self):
        pass

class Parallel:

    MAX_PROCESSES = 78
    MAX_QUEUE = MAX_PROCESSES+50

    def __init__(self):
        # empty job queue
        self._queue = JoinableQueue(Parallel.MAX_QUEUE)
        self.workers: List[Parallel.Consumer] = []
        self.shared_session: Optional[Session] = None

    def import_virus_sample(self, session: Session, sample: VirusSample):
        # do this synchronously
        experiment = vcm.create_or_get_experiment(session, sample)
        host_sample = vcm.create_or_get_host_sample(session, sample)
        sequencing_project = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_or_get_sequence(session, sample, virus_id, experiment, host_sample,
                                              sequencing_project)

        if not self.workers:
            self.shared_session = database_tom.get_session()
            # prepare workers
            for _ in range(Parallel.number_of_processes()):
                worker = Parallel.Consumer(self._queue, virus.nucleotide_variant_aligner(), self.shared_session)
                self.workers.append(worker)
                worker.start()

        # schedule nucleotide variants to be called asynchronously
        sample.on_before_multiprocessing()
        self._queue.put([sample, sequence.sequence_id])
        logger.info(f'nucleotide variant calling for sequence {sample.internal_id()} scheduled. Queue size: {self._queue.qsize()}')

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


def import_virus(session: Session, virus: VirusSource):
    return vcm.create_or_get_virus(session, virus).virus_id


def try_import_virus_sample(sample: VirusSample):
    global successful_imports
    try:
        database_tom.try_py_function(import_method.import_virus_sample, sample)
        successful_imports += 1
    except:
        logger.exception(f'exception occurred while working on virus sample {sample}')


virus = COGUKSarsCov2()
# IMPORT VIRUS TAXON DATA
virus_id = database_tom.try_py_function(import_virus, virus)

# # IMPORT VIRUS SEQUENCES
logger.info(f'importing virus sequences and related tables')
import_method = Parallel()
successful_imports = 0

for s in virus.virus_samples():
    try_import_virus_sample(s)

logger.info('main process completed')
import_method.tear_down()
logger.info(f'successful imports: {successful_imports} (not reliable when parallel processing)')
