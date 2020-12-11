import os
import pickle
from loguru import logger
from typing import Optional, List
from queuable_tasks import max_number_of_workers
from vcm import vcm as vcm
from data_sources.coguk_sars_cov_2.sample import COGUKSarsCov2Sample
from data_sources.coguk_sars_cov_2.virus import COGUKSarsCov2
from multiprocessing import JoinableQueue, cpu_count, Process
from sqlalchemy.orm.session import Session
from db_config import read_db_import_configuration as import_config, database_tom
from data_sources.ncbi_any_virus.settings import known_settings as ncbi_known_settings
from locations import get_local_folder_for, FileType, remove_file
from pipeline_nuc_variants__annotations__aa import sequence_aligner
import stats_module

sc2_chromosome = ncbi_known_settings["sars_cov_2"]["chromosome_name"]
sc2_annotations_file_path = ncbi_known_settings["sars_cov_2"]["annotation_file_path"]
sc2_snpeff_db_name = ncbi_known_settings["sars_cov_2"]["snpeff_db_name"]


reference_sequence = None
virus: Optional[COGUKSarsCov2] = None
virus_id: Optional[int] = None
import_method = None
successful_imports: Optional[int] = None


def main_pipeline_part_3(session: database_tom.Session, sample, db_sequence_id):
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
        sequence = vcm.create_and_get_sequence(session, sample, virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)
        main_pipeline_part_3(session, sample, sequence.sequence_id)

    def tear_down(self):
        pass

class Parallel:

    MAX_PROCESSES = max_number_of_workers(50)

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
        sequence = vcm.create_and_get_sequence(session, sample, virus_id, experiment_id, host_sample_id,
                                               sequencing_project_id)

        if not self.workers:
            global reference_sequence
            reference_sequence = virus.reference_sequence()
            # prepare workers
            for _ in range(Parallel.number_of_processes()):
                worker = Parallel.Consumer(self._queue, virus.reference_sequence(), database_tom.get_session())
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
                            database_tom.rollback(self.shared_session)
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


def import_virus(session: Session, _virus: COGUKSarsCov2):
    return vcm.create_or_get_virus(session, _virus)


def run(from_sample: Optional[int] = None, to_sample: Optional[int] = None):
    global virus, virus_id, import_method, successful_imports
    db_params: dict = import_config.get_database_config_params()
    database_tom.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])
    virus = COGUKSarsCov2()
    # IMPORT VIRUS TAXON DATA
    virus_id = database_tom.try_py_function(import_virus, virus)

    # update last import date
    database_tom.try_py_function(vcm.update_db_metadata, virus_id, 'COG-UK')

    # # IMPORT VIRUS SEQUENCES
    logger.info(f'importing virus sequences and related tables')
    import_method = Parallel()
    successful_imports = 0

    def try_import_virus_sample(sample: COGUKSarsCov2Sample):
        global successful_imports
        try:
            database_tom.try_py_function(import_method.import_virus_sample, sample)
            successful_imports += 1
        except:
            logger.exception(f'exception occurred while working on virus sample {sample.primary_accession_number()}')

    # total_s = 2
    for s in virus.virus_samples(virus_id, from_sample, to_sample):
        if not s.nucleotide_sequence():
            logger.info(f'sample {s.primary_accession_number()} skipped because nucleotide sequence is empty or null')
            continue
        # if total_s > 0:
        #     total_s -= 1
        # else:
        #     break
        try_import_virus_sample(s)

    logger.info('main process completed')
    import_method.tear_down()
    logger.info(f'successful imports: {successful_imports} (not reliable when parallel processing)')
