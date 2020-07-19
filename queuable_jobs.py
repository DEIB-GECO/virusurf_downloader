from multiprocessing import JoinableQueue, cpu_count, Process
import os
from typing import List, Optional
from loguru import logger


class Job:

    def execute(self):
        raise NotImplementedError()


class Worker(Process):
    def __init__(self, jobs: JoinableQueue):
        super().__init__()
        self.jobs = jobs
        logger.info('worker started')

    def get_a_job(self) -> Optional[Job]:
        return self.jobs.get(block=True, timeout=None)

    def run(self):
        while True:
            job = self.get_a_job()
            if job is None:
                logger.info('shutting down a consumer')
                self.jobs.task_done()
                self.release_resources()
                break
            else:
                try:
                    job.execute()
                except Exception:
                    logger.exception('An exception reached the Worker loop. Worker is being disposed.')
                    self.release_resources()
                    break
                finally:
                    self.jobs.task_done()

    def release_resources(self):
        self.jobs = None


class Boss:
    # noinspection PyPep8Naming
    def __init__(self, jobs_queue_capacity: int, workers_num: int, WorkerClass: Worker.__class__ = Worker):
        # empty job queue
        self._queue = JoinableQueue(maxsize=jobs_queue_capacity)
        logger.info(f'Queue size set to accept at most {jobs_queue_capacity} before pausing job assignment.')
        self.WorkerClass = WorkerClass
        self.workers_num = workers_num

    _workers = []

    def wake_up_workers(self):
        self._workers: List[Worker] = [self.WorkerClass(self._queue) for _ in range(self.workers_num)]
        for worker in self._workers:
            worker.start()

    def assign_job(self, job: Job):
        self._queue.put(job)

    def stop_workers(self):
        for _ in self._workers:
            self._queue.put(None, block=True, timeout=None)
        logger.info('waiting all workers to finish')
        self._queue.join()
        logger.info('all processes_finished')


def max_number_of_workers(user_defined_max_processes):
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
    max_p.append(user_defined_max_processes)
    used_processes = min(max_p)
    if used_processes < user_defined_max_processes:
        logger.warning(f'maximum number of CPUs or usable process reached. At most {used_processes} will be used.')
    return used_processes
