import os
from typing import Optional, Callable, List

from loguru import logger

import database_tom
import stats_module
from data_sources.virus_sample import VirusSample
from database_tom import Session
import vcm
from data_sources.gisaid_sars_cov_2.virus import GISAIDSarsCov2
from multiprocessing import JoinableQueue, cpu_count, Process


class Sequential:

    def __init__(self, virus_id: int):
        self.virus_id = virus_id
        self.aligner: Optional[Callable] = None
        self.reference_sample: Optional[VirusSample] = None

    def import_virus_sample(self, session: Session, sample: VirusSample):
        experiment_id = vcm.create_or_get_experiment(session, sample)
        host_sample_id = vcm.create_or_get_host_sample(session, sample)
        sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
        sequence = vcm.create_and_get_sequence(session, sample, self.virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)
        vcm.create_annotation_and_aa_variants(session, sample, sequence, self.reference_sample)
        stats_module.completed_sample()

    def tear_down(self):
        pass


# class Parallel:
#
#     MAX_PROCESSES = 3
#
#     def __init__(self, virus_id: int):
#         self.virus_id = virus_id
#         self.reference_sample: Optional[VirusSample] = None
#         # empty job queue
#         queue_size = self.number_of_processes()+10
#         self._queue = JoinableQueue(maxsize=queue_size)
#         logger.info(f'Queue size set to accept at most {queue_size} before pausing the producer.')
#         self.workers: List[Parallel.Consumer] = []
#         self.shared_session: Optional[Session] = None
#
#     def import_virus_sample(self, session: Session, sample):
#         # do this synchronously
#         experiment_id = vcm.create_or_get_experiment(session, sample)
#         host_sample_id = vcm.create_or_get_host_sample(session, sample)
#         sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
#         sequence = vcm.create_and_get_sequence(session, sample, self.virus_id, experiment_id, host_sample_id,
#                                               sequencing_project_id)
#         vcm.create_annotation_and_aa_variants(session, sample, sequence, self.reference_sample)
#
#         if sample.is_reference():
#             self.reference_sample = sample
#             self.shared_session = database_tom.get_session()
#             # prepare workers
#             for _ in range(Parallel.number_of_processes()):
#                 worker = Parallel.Consumer(self._queue, sample.nucleotide_var_aligner(), self.shared_session)
#                 self.workers.append(worker)
#                 worker.start()
#             logger.info('reference sequence imported')
#         else:
#             # schedule nucleotide variants to be called asynchronously
#             sample.on_before_multiprocessing()
#             self._queue.put([sample, sequence.sequence_id])
#             logger.info(f'nucleotide variant calling for sequence {sample.internal_id()} scheduled')
#
#     class Consumer(Process):
#         def __init__(self, jobs: JoinableQueue, aligner: Callable, shared_session: Session):
#             super().__init__()
#             self.jobs = jobs
#             self.aligner: Optional[Callable] = aligner
#             self.shared_session: Session = shared_session
#             logger.info('worker started')
#
#         def run(self):
#             while True:
#                 job = self.jobs.get(block=True, timeout=None)
#                 if job is None:
#                     logger.info('shutting down a consumer')
#                     self.jobs.task_done()
#                     break
#                 else:
#                     sample, sequence_id = job
#                     try:
#                         vcm.create_nucleotide_variants_and_impacts(self.shared_session, sample, sequence_id, self.aligner)
#                         self.shared_session.commit()
#                         logger.info(f'nucleotides of sequence {sample.internal_id()} imported')
#                         stats_module.completed_sample()
#                     except:
#                         logger.exception(f'unknown exception while calling nucleotides on sample {sample.internal_id()}')
#                         database_tom.rollback(self.shared_session)
#                     self.jobs.task_done()
#
#     @staticmethod
#     def number_of_processes():
#         max_p = []
#         # get max allowed processes
#         try:
#             max_p.append(len(os.sched_getaffinity(0)))
#         except AttributeError:
#             pass    # method not available on this system (can happen on Windows and macOS)
#         # get number of CPUs
#         try:
#             max_p.append(cpu_count())
#         except NotImplementedError:
#             pass
#         # take the smallest number of above two and Parallel.MAX_PROCESSES
#         max_p.append(Parallel.MAX_PROCESSES)
#         used_processes = min(max_p)
#         if used_processes < Parallel.MAX_PROCESSES:
#             logger.warning(f'maximum number of CPUs or usable process reached. At most {used_processes} will be used.')
#         return used_processes
#
#     def tear_down(self):
#         # terminate workers
#         for _ in self.workers:
#             self._queue.put(None, block=True, timeout=None)
#         logger.info('waiting all workers to finish')
#         self._queue.join()
#         logger.info('all processes_finished')
#         self.shared_session.close()


def run(from_sample: Optional[int] = None, to_sample: Optional[int] = None):

    def import_virus(session: Session, virus: GISAIDSarsCov2):
        return vcm.create_or_get_virus(session, virus)

    def try_import_virus_sample(sample: VirusSample):
        nonlocal successful_imports, import_method
        try:
            database_tom.try_py_function(import_method.import_virus_sample, sample)
            successful_imports += 1
        except:
            logger.exception(f'exception occurred while working on virus sample {sample.internal_id()}')

    # IMPORT VIRUS TAXON DATA
    virus = GISAIDSarsCov2()
    virus_id = database_tom.try_py_function(import_virus, virus)

    # # IMPORT VIRUS SEQUENCES
    logger.info(f'importing virus sequences and related tables')
    import_method = Sequential(virus_id)
    successful_imports = 0

    for s in virus.virus_samples(virus_id, from_sample, to_sample):
        try_import_virus_sample(s)

    logger.info('main process completed')
    import_method.tear_down()
    logger.info(f'successful imports: {successful_imports} (not reliable when parallel processing)')
