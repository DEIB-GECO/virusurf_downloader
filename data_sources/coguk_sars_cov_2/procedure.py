import os
import pickle
from time import sleep
from loguru import logger
from typing import Optional, List
from queuable_tasks import max_number_of_workers
from vcm import vcm as vcm
from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from data_sources.coguk_sars_cov_2.virus import COGUKSarsCov2
from multiprocessing import JoinableQueue, cpu_count, Process
from sqlalchemy.orm.session import Session
from db_config import read_db_import_configuration as import_config, database
from data_sources.ncbi_any_virus.settings import known_settings as ncbi_known_settings
from locations import get_local_folder_for, FileType, remove_file
from nuc_aa_pipeline import sequence_aligner
import stats_module

sc2_chromosome = ncbi_known_settings["sars_cov_2"]["chromosome_name"]
sc2_annotations_file_path = ncbi_known_settings["sars_cov_2"]["annotation_file_path"]
sc2_snpeff_db_name = ncbi_known_settings["sars_cov_2"]["snpeff_db_name"]


reference_sequence = None
virus: Optional[COGUKSarsCov2] = None
virus_id: Optional[int] = None
import_method = None


def main_pipeline_part_3(session: database.Session, sample, db_sequence_id):
    file_path = get_local_folder_for(virus.name, FileType.Annotations) + str(sample.internal_id()) + ".pickle"
    try:
        if not os.path.exists(file_path):
            annotations_and_nuc_variants = sequence_aligner(
                sample.internal_id(),
                reference_sequence,
                sample.nucleotide_sequence(),
                sc2_chromosome,
                sc2_annotations_file_path,
                sc2_snpeff_db_name)
            with open(file_path, mode='wb') as cache_file:
                pickle.dump(annotations_and_nuc_variants, cache_file, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(file_path, mode='rb') as cache_file:
                annotations_and_nuc_variants = pickle.load(cache_file)
        annotations, nuc_variants = annotations_and_nuc_variants
        for ann in annotations:
            vcm.create_annotation_and_amino_acid_variants(session, db_sequence_id, *ann)
        for nuc in nuc_variants:
            vcm.create_nuc_variants_and_impacts(session, db_sequence_id, nuc)
        stats_module.completed_sample(sample.primary_accession_number())
    except Exception as e:
        if str(e).endswith("sequence contains letters not in the alphabet"):
            logger.warning(f"sample {sample.primary_accession_number()} skipped because sequence contains letter not in "
                           f"the alphabet")
        else:
            logger.exception(f'exception occurred during pipeline_part_3 of sample {sample.primary_accession_number()}. '
                         f'Doing rollback of insertion of variants and annotations + deletion of cache')
        remove_file(file_path)
        raise e


class Sequential:

    def __init__(self):
        global reference_sequence
        reference_sequence = virus.reference_sequence()

    def import_virus_sample(self, session: Session, sample: COGUKSarsCov2Sample):
        experiment_id = vcm.create_or_get_experiment(session, sample)
        host_specie_id = vcm.create_or_get_host_specie(session, sample)
        host_sample_id = vcm.create_or_get_host_sample(session, sample, host_specie_id)
        sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
        sequence, nucleotide_seq = vcm.create_and_get_sequence(session, sample, virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)
        main_pipeline_part_3(session, sample, sequence.sequence_id)

    def tear_down(self):
        pass


class Parallel:

    MAX_PROCESSES = max_number_of_workers(12)

    def __init__(self):
        # empty job queue
        queue_size = self.number_of_processes()+10
        self._queue = JoinableQueue(queue_size)
        logger.info(f'queue set to accapt at most {queue_size} jobs before pausing the producer')
        self.workers: List[Parallel.Consumer] = []

    def import_virus_sample(self, session: Session, sample: COGUKSarsCov2Sample):
        # do this synchronously
        experiment_id = vcm.create_or_get_experiment(session, sample)
        host_specie_id = vcm.create_or_get_host_specie(session, sample)
        host_sample_id = vcm.create_or_get_host_sample(session, sample, host_specie_id)
        sequencing_project_id = vcm.create_or_get_sequencing_project(session, sample)
        sequence, nucleotide_seq = vcm.create_and_get_sequence(session, sample, virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)

        if not self.workers:
            global reference_sequence
            reference_sequence = virus.reference_sequence()
            # prepare workers
            for _ in range(Parallel.number_of_processes()):
                worker = Parallel.Consumer(self._queue, virus.reference_sequence(), database.get_session())
                self.workers.append(worker)
                worker.start()

        # schedule nucleotide variants to be called asynchronously
        sample.on_before_multiprocessing()
        self._queue.put([sample, sequence.sequence_id])
        # logger.debug(f'nuc_var for sequence {sample.primary_accession_number()} scheduled\tQueue size: {self._queue.qsize()}')

    class Consumer(Process):
        def __init__(self, jobs: JoinableQueue, refseq: str, shared_session: Session):
            super().__init__()
            self.jobs = jobs
            self.refseq = refseq
            self.shared_session: Session = shared_session
            logger.info('worker started')

        def run(self):
            try:
                while True:
                    job = self.jobs.get(block=True, timeout=None)
                    if job is None:
                        logger.info('shutting down a consumer')
                        self.jobs.task_done()
                        break
                    else:
                        sample, sequence_id = job
                        try:
                            main_pipeline_part_3(self.shared_session, sample, sequence_id)
                            self.shared_session.commit()
                        except:
                            database.rollback(self.shared_session)
                        finally:
                            self.jobs.task_done()
            finally:
                self.shared_session.close()
                self.jobs = None

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


def run(from_sample: Optional[int] = None, to_sample: Optional[int] = None):
    global virus, virus_id, import_method
    db_params: dict = import_config.get_database_config_params()
    database.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])
    virus = COGUKSarsCov2()
    # IMPORT VIRUS TAXON DATA
    virus_id = database.try_py_function(vcm.create_or_get_virus, virus)

    # update last import date
    database.try_py_function(vcm.update_db_metadata, virus_id, 'COG-UK')

    # find outdated and new samples from source (some sequences can be updated, so the sets are not necessarily disjoint)
    logger.warning(
        "Current implementation of deltas for COG-UK uses more than 10 GB of RAM to cache query results and save time.\n"
        "IF YOUR SYSTEM CAN'T PROVIDE MORE THAN 10 GB OF RAM, STOP THE PROCESS NOW.\n"
        "The program will resume in 15 seconds")
    try:
        sleep(15)
    except KeyboardInterrupt:
        return
    id_outdated_sequences, id_new_sequences = virus.deltas()
    logger.warning('Check deltas.. The program will resume in 30 seconds.')
    try:
        sleep(30)
    except KeyboardInterrupt:
        return

    # select range
    if from_sample is not None and to_sample is not None:
        id_new_sequences = {id_new_sequences.pop() for i in range(to_sample-from_sample)}

    # initialize statistics module
    stats_module.schedule_samples(
        stats_module.StatsBasedOnIds(id_new_sequences, True, virus_id, ['COG-UK']))

    # remove outdated sequences
    logger.info(f'removing outdated sequences')
    database.try_py_function(vcm.remove_sequence_and_meta_list, primary_sequence_accession_id=id_outdated_sequences)
    stats_module.removed_samples(id_outdated_sequences)
    for _id in id_outdated_sequences:
        file_path = get_local_folder_for(virus.name, FileType.Annotations) + str(_id).replace('/', '-') + ".pickle"
        remove_file(file_path)

    # prepare multiprocessing
    logger.info(f'importing virus sequences and related tables')
    import_method = Parallel()

    vcm.DBCache.commit_changes()
    for s in virus.get_sequences_of_updated_source(filter_accession_ids=id_new_sequences):
        if not s.nucleotide_sequence():
            logger.info(f'sample {s.primary_accession_number()} skipped because nucleotide sequence is empty or null')
            continue
        try:
            database.try_py_function(import_method.import_virus_sample, s)
            vcm.DBCache.commit_changes()
        except:
            logger.exception(f'exception occurred while working on virus sample {s.primary_accession_number()}')
            vcm.DBCache.rollback_changes()

    logger.info('main process completed')
    import_method.tear_down()

    # remove leftovers of failed samples
    try:
        metadata_samples_to_remove: set = stats_module.get_scheduled_not_completed()
        if len(metadata_samples_to_remove) > 0:
            logger.info(
                f"Removing metadata leftovers of imports that failed during variant/annotation calling or metadata"
                f" ({len(metadata_samples_to_remove)} samples)")

            metadata_samples_to_remove_as_string: list = [str(x) for x in metadata_samples_to_remove]
            logger.trace("Accession id of failed imports:\n"
                         f"{metadata_samples_to_remove_as_string}")
            logger.info("Deleting leftovers in database")
            database.try_py_function(vcm.remove_sequence_and_meta_list,
                                     primary_sequence_accession_id=metadata_samples_to_remove_as_string)
    except:
        logger.exception(
            "Removal of metadata leftovers in the DB of the samples that failed was not successful.")

