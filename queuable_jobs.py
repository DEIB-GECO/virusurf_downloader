from multiprocessing import JoinableQueue, cpu_count, Process
import os
from queue import Empty
from typing import List, Optional
from loguru import logger


class Job:

    def execute(self):
        """
        It's advisable to handle KeyboardInterrupt exceptions by doing any necessary clean-up operation before returning.
        If there are no clean-up operations, you may want to pause the process for one second to let the main thread have
        the time to cancel scheduled jobs before eager workers get them.
        """
        raise NotImplementedError()


class Worker(Process):
    def __init__(self, jobs: JoinableQueue):
        super().__init__()
        self.jobs = jobs
        logger.info('worker started')

    def get_a_job(self) -> Optional[Job]:
        """
        :return: a Job instance. None objects are reserved for use by the Boss.
        """
        return self.jobs.get(block=True, timeout=None)

    def run(self):
        """
        KeyboardInterrupt are handled only when occurring in the two most common conditions:
        - A worker is waiting to receive a job.
        - A worker is executing a supposedly long task and it's interrupted before it is completed.
        In both cases, the queue is maintained in a consistent state and the worker process won't hang. The best
        way to terminate the processes is by catching the KeyboardInterrupt in the main thread and use the prepared
        methods to discard the remaining jobs and stop the workers.
        """
        while True:
            try:
                job = self.get_a_job()
            except KeyboardInterrupt:
                continue   # Ctrl+C arrived while waiting to receive a job
            else:
                if job is None:
                    self.jobs.task_done()
                    self.release_resources()
                    break
                else:
                    try:
                        job.execute()
                        self.jobs.task_done()
                    except KeyboardInterrupt:
                        self.jobs.task_done()
                        logger.info('KeyboardInterrupt caused a worker to interrupt a job before completing it. '
                                    'You may want to handle this exception inside the Job definition instead.')
                    except Exception:
                        logger.exception('An exception reached the Worker loop. Worker is being disposed.')
                        self.jobs.task_done()
                        self.release_resources()
                        break
        logger.info('worker left')

    def release_resources(self):
        self.jobs = None


class Boss:
    # noinspection PyPep8Naming
    def __init__(self, jobs_queue_capacity: int, workers_num: int, WorkerClass: Worker.__class__ = Worker):
        # empty job queue
        self._queue = JoinableQueue(maxsize=jobs_queue_capacity)
        logger.info(f'Queue size set to accept at most {jobs_queue_capacity} before pausing job assignment.')
        self.WorkerClass = WorkerClass
        self.workers_num = max_number_of_workers(workers_num)

    _workers = []

    def wake_up_workers(self):
        self._workers: List[Worker] = [self.WorkerClass(self._queue) for _ in range(self.workers_num)]
        for worker in self._workers:
            worker.start()

    def assign_job(self, job: Job):
        self._queue.put(job)

    def stop_workers(self):
        logger.info('waiting all workers to finish')
        for _ in self._workers:
            self._queue.put(None, block=True, timeout=None)
        self._queue.join()
        logger.info('all processes_finished')

    def discard_left_jobs(self):
        while not self._queue.empty():
            try:
                self._queue.get(False)
            except Empty:
                continue
            self._queue.task_done()

    def number_of_waiting_jobs(self):
        return self._queue.qsize()


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
