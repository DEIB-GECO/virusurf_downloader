from datetime import datetime
from typing import Optional
from loguru import logger
import stats_module
from data_sources.virus_sample import VirusSample
from db_config.database import Session
from vcm import vcm
from data_sources.gisaid_sars_cov_2.virus import GISAIDSarsCov2
from time import sleep
from tqdm import tqdm
from db_config import read_db_import_configuration as import_config, database
from logger_settings import send_message


class Sequential:

    def __init__(self, virus_id: int):
        self.virus_id = virus_id

    def import_virus_sample(self, session: Session, sample: VirusSample):
        experiment_id = vcm.create_or_get_experiment(session, sample)
        host_specie_id = vcm.create_or_get_host_specie(session, sample)
        host_sample_id = vcm.create_or_get_host_sample(session, sample, host_specie_id)
        sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
        sequence, nucleotide_seq = vcm.create_and_get_sequence(session, sample, self.virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)
        vcm.create_annotation_and_aa_variants(session, sample, sequence, None)
        stats_module.completed_sample(sample.primary_accession_number())

    def update_virus_sample(self, session: Session, sample: VirusSample, changes: dict):
        host_sample_id = None
        sequencing_project_id = None
        sequence = None
        if changes["host_sample"] is True:
            host_specie_id = vcm.create_or_get_host_specie(session, sample)
            host_sample_id = vcm.create_or_get_host_sample(session, sample, host_specie_id)
        if changes["sequencing_project"] is True:
            sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
        if changes["sequence"] is True or changes["host_sample"] is True or changes["sequencing_project"] is True:
            sequence = vcm.update_and_get_sequence(session, sample,
                                                   experiment_id=None,
                                                   host_sample_id=host_sample_id,
                                                   sequencing_project_id=sequencing_project_id)
        if changes["annotations"] is True:
            if sequence is None:
                sequence = vcm.get_sequence(session, sample, self.virus_id) # need the sequence_id
            vcm.remove_annoations_and_aa_variants(session=session, sequence_id=sequence.sequence_id)
            vcm.create_annotation_and_aa_variants(session, sample, sequence, None)

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
#             self.shared_session = database.get_session()
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
#                         database.rollback(self.shared_session)
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
    db_params: dict = import_config.get_database_config_params()
    database.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])

    # IMPORT VIRUS TAXON DATA
    virus = GISAIDSarsCov2()
    virus_id = database.try_py_function(vcm.create_or_get_virus, virus)
    database.try_py_function(vcm.update_db_metadata, virus_id, 'GISAID')

    # COMPUTE DELTAS
    acc_ids_sequences_to_remove, acc_id_sequences_to_import, sequences_to_update = virus.deltas()
    acc_ids_sequences_to_update = sequences_to_update.keys()
    logger.warning('Check deltas.. The program will resume in 30 seconds.')
    sleep(30)

    # MIND FROM_SAMPLE/TO_SAMPLE
    if from_sample is not None and to_sample is not None:
        count_el_to_import = abs(to_sample - from_sample)
        count_el_to_import = min(count_el_to_import, len(acc_id_sequences_to_import))
        count_el_to_ignore = len(acc_id_sequences_to_import) - count_el_to_import
        try:
            for i in range(count_el_to_ignore):
                acc_id_sequences_to_import.pop()
        except KeyError:
            pass    # if pop on empty set

    # create pipeline_event (will be inserted later)
    pipeline_event = database.PipelineEvent(
        event_date=datetime.now().strftime("%Y-%m-%d"),
        event_name=f'GISAID sars_cov_2 sequences update',
        removed_items=len(acc_ids_sequences_to_remove),
        changed_items=len(acc_ids_sequences_to_update),     # provisional - (sqlalchemy wants a value at obj creation)
        added_items=len(acc_id_sequences_to_import)         # provisional - (sqlalchemy wants a value at obj creation)
    )
    changed_items = 0
    added_items = 0

    stats_module.schedule_samples(stats_module.StatsBasedOnIds(acc_id_sequences_to_import, True, virus_id, ['GISAID']))

    logger.info('Removing outdated sequences...')
    # REMOVE OUTDATED SEQUENCES
    database.try_py_function(vcm.remove_sequence_and_meta_list, acc_ids_sequences_to_remove, None)
    stats_module.removed_samples(acc_ids_sequences_to_remove)

    # IMPORT NEW/CHANGED SEQUENCES
    vcm.DBCache.commit_changes()
    logger.info(f'Importing virus sequences and related tables...')
    import_method = Sequential(virus_id)

    progress = tqdm(total=len(acc_id_sequences_to_import)+len(acc_ids_sequences_to_update))
    for sample in virus.get_sequences_of_updated_source():
        sample_accession_id = sample.primary_accession_number()
        try:
            if sample_accession_id in acc_id_sequences_to_import:
                # import sample from scratch
                database.try_py_function(import_method.import_virus_sample, sample)
                vcm.DBCache.commit_changes()
                progress.update()
                added_items += 1
            elif sample_accession_id in acc_ids_sequences_to_update:
                # update values inside the database
                changes_in_sequence = sequences_to_update[sample_accession_id]
                database.try_py_function(import_method.update_virus_sample, sample, changes_in_sequence)
                vcm.DBCache.commit_changes()
                progress.update()
                changed_items += 1
        except KeyboardInterrupt:
            logger.info("main loop interrupted by the user")
            break
        except:
            logger.exception(f'exception occurred while working on virus sample {sample.internal_id()}')
            vcm.DBCache.rollback_changes()

    logger.info('main loop completed')
    import_method.tear_down()

    logger.info('Removal of unused database objects...')
    database.try_py_function(vcm.clean_objects_unreachable_from_sequences)

    if len(acc_id_sequences_to_import) - added_items > 100:
        send_message(f"GISAID importer can have a bug. import of {len(acc_id_sequences_to_import) - added_items} out of"
                     f" {len(acc_id_sequences_to_import)} failed.")

    pipeline_event.changed_items = changed_items
    pipeline_event.added_items = added_items
    database.try_py_function(vcm.insert_data_update_pipeline_event, pipeline_event)
