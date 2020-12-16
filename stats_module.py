import sys
from datetime import datetime
from queue import Empty
from typing import Collection, List, Optional
from multiprocessing import Queue, Value
from loguru import logger
from db_config import database_tom
from vcm import vcm
from abc import ABC
from collections import Counter


class StatsType(ABC):
    def __init__(self):
        self._start_time = datetime.now()
        self._is_primary_acc_id = True
        self._virus_db_id = None
        self._sources = None

    def _measure_delta_in_db(self, virus_db_id: int, is_primary_acc_id: bool, sources: Optional[List[str]] = None):
        self._is_primary_acc_id = is_primary_acc_id
        self._virus_db_id = virus_db_id
        self._sources = sources
        if is_primary_acc_id:
            _sequences_in_db_at_start = database_tom.try_py_function(vcm.sequence_primary_accession_ids, virus_db_id,
                                                                     sources)
        else:
            _sequences_in_db_at_start = database_tom.try_py_function(vcm.sequence_alternative_accession_ids,
                                                                     virus_db_id, sources)
        self._sequences_in_db_at_start = Counter(_sequences_in_db_at_start)

    def _performance_message(self, num_processed_samples: int):
        # performance
        end_time = datetime.now()
        elapsed_time = end_time - self._start_time
        add_to_message = f"Elapsed Time {elapsed_time}\n"
        speed = elapsed_time / num_processed_samples if num_processed_samples != 0 else "--"
        add_to_message += f"Speed {speed} sec./samples\n"
        return add_to_message

    def check_samples_imported(self, *args, **kwargs):
        raise NotImplementedError

    def completed_sample(self, *args, **kwargs):
        raise NotImplementedError

    def removed_samples(self, *args, **kwargs):
        raise NotImplementedError


class StatsBasedOnTotals(StatsType):
    def __init__(self, number_of_scheduled_samples: int, virus_db_id: Optional[int] = None, sources: Optional[List[str]] = None):
        super().__init__()
        self._number_scheduled_samples = number_of_scheduled_samples
        self._number_completed_samples = Value('i', 0)  # multiprocess-safe int initialized to 0
        self._number_removed_samples = Value('i', 0)  # multiprocess-safe int initialized to 0
        message = f"\n" \
                  f"STATS MODULE:\n" \
                  f"scheduled import of {number_of_scheduled_samples} samples"
        logger.info(message)
        if virus_db_id is not None:
            self._measure_delta_in_db(virus_db_id, True, sources)

    def completed_sample(self):
        with self._number_completed_samples.get_lock():
            self._number_completed_samples.value += 1

    def removed_samples(self, number):
        with self._number_removed_samples.get_lock():
            self._number_removed_samples.value += number

    def check_samples_imported(self):
        # overwrite Value object with int for a simpler handling
        self._number_completed_samples = self._number_completed_samples.value
        self._number_removed_samples = self._number_removed_samples.value

        message = f"\n" \
                  f"STATS MODULE:\n" \
                  f"Scheduled import of {self._number_scheduled_samples} samples\n" \
                  f"Removed {self._number_removed_samples} samples\n" \
                  f"Completed {self._number_completed_samples} samples\n"
        warn = False

        # find errors
        scheduled_not_completed = self._number_scheduled_samples - self._number_completed_samples
        completed_not_scheduled = self._number_completed_samples - self._number_scheduled_samples
        if scheduled_not_completed > 0:
            warn = True
            message += f'Failed (scheduled but not completed) {scheduled_not_completed} samples\n'
        if completed_not_scheduled > 0:
            warn = True
            message += f'Wrongly processed (completed but not scheduled) {completed_not_scheduled} samples\n'

        # check against DB
        if self._virus_db_id is not None:
            current_sequences_in_db = database_tom.try_py_function(
                vcm.sequence_primary_accession_ids, self._virus_db_id, self._sources)
            current_sequences_in_db = Counter(current_sequences_in_db)
            inserted_in_db = current_sequences_in_db - self._sequences_in_db_at_start
            removed_from_db = self._sequences_in_db_at_start - current_sequences_in_db
            message += f"Sequences in DB:\n" \
                       f"\tprevious number {sum(self._sequences_in_db_at_start.values())}.\n" \
                       f"\tcurrent number {sum(current_sequences_in_db.values())}\n" \
                       f"\tdifference: {sum(inserted_in_db.values())} new - {sum(removed_from_db.values())} missing from source and deleted = {sum(inserted_in_db.values()) - sum(removed_from_db.values())}\n"
            # check duplicates
            num_duplcates_at_start = sum(self._sequences_in_db_at_start.values()) - len(set(self._sequences_in_db_at_start))
            num_current_duplicates = sum(current_sequences_in_db.values()) - len(set(current_sequences_in_db))
            if num_duplcates_at_start > 0 or num_current_duplicates > 0:
                warn = True
                message += f"Duplicated accession_ids in DB:\n" \
                           f"\tprevious number of duplicates: {num_duplcates_at_start}\n" \
                           f"\tcurrent number of duplicates: {num_current_duplicates}\n" \
                           f"\tdetail of current duplicated accession_ids: {sorted(list(current_sequences_in_db - Counter(set(current_sequences_in_db))))}\n"

            num_completed_not_inserted = self._number_completed_samples - sum(inserted_in_db.values()) - self._number_removed_samples
            if num_completed_not_inserted > 0:
                warn = True
                message += f'Number of samples completed and not imported (or imported twice) into the DB: {num_completed_not_inserted}\n'
            elif num_completed_not_inserted < 0:
                warn = True
                message += f'Number of samples imported into the DB but incomplete: {-num_completed_not_inserted}\n'

            # performance
            message += self._performance_message(self._number_completed_samples)

        if not warn:
            logger.info(message)
        else:
            logger.warning(message)


class StatsBasedOnIds(StatsType):
    def __init__(self,samples_id: Collection[str], is_primary_acc_id: bool, virus_db_id: Optional[int] = None, sources: Optional[List[str]] = None):
        super().__init__()
        self._scheduled_samples_acc_id = [str(x) for x in samples_id]
        message = f"\n" \
                  f"STATS MODULE:\n" \
                  f"scheduled import of {len(self._scheduled_samples_acc_id)} samples"
        if virus_db_id is not None:
            self._measure_delta_in_db(virus_db_id, is_primary_acc_id, sources)
        if len(set(samples_id)) == len(samples_id):
            logger.info(message)
        else:
            set_scheduled = set(self._scheduled_samples_acc_id)
            _duplicated_accession_id = list(self._scheduled_samples_acc_id)
            for x in set_scheduled:
                _duplicated_accession_id.remove(x)
            difference = len(_duplicated_accession_id)
            message = f'\n{difference} of the scheduled samples appear to have '
            message += 'primary ' if is_primary_acc_id else 'alternative '
            message += " duplicated accession ids. Ids are\n"
            message += f"{sorted(_duplicated_accession_id)}"
            logger.error(message)
            sys.exit(1)
        self._completed_samples_acc_id = Queue()
        self._removed_samples_acc_id = Queue()

    def completed_sample(self, sample_acc_id: str):
        self._completed_samples_acc_id.put_nowait(sample_acc_id)

    def removed_samples(self, samples_acc_id: List[str]):
        for i in samples_acc_id:
            self._removed_samples_acc_id.put_nowait(i)

    def check_samples_imported(self):
        if self._completed_samples_acc_id is None:
            logger.error('STATS MODULE: check_samples_imported called before add_samples. Stats cannot be produced.')
        else:
            warn = False
            message = f"\n" \
                      f"STATS MODULE:\n"
            _scheduled_samples_acc_id = set(self._scheduled_samples_acc_id)
            message += f"Scheduled import of {len(_scheduled_samples_acc_id)} samples\n"
            completed_samples_acc_id = self._queue_to_set(self._completed_samples_acc_id)
            removed_samples_acc_id = self._queue_to_set(self._removed_samples_acc_id)
            self._completed_samples_acc_id.close()
            self._removed_samples_acc_id.close()
            message += f"Removed {len(removed_samples_acc_id)} samples\n"
            message += f"Completed {len(completed_samples_acc_id)} samples\n"

            # find errors
            scheduled_not_completed = _scheduled_samples_acc_id - completed_samples_acc_id
            completed_not_scheduled = completed_samples_acc_id - _scheduled_samples_acc_id
            if len(scheduled_not_completed) > 0:
                warn = True
                message += f'Failed (scheduled but not completed) {len(scheduled_not_completed)} samples\n' \
                           f'\tsamples id: {scheduled_not_completed}\n'
            if len(completed_not_scheduled) > 0:
                warn = True
                message += f'Wrongly processed (completed but not scheduled) {len(completed_not_scheduled)} samples\n' \
                           f'\tsamples id: {completed_not_scheduled}\n'

            # check against DB
            if self._virus_db_id is not None:
                if self._is_primary_acc_id:
                    current_sequences_in_db = database_tom.try_py_function(
                        vcm.sequence_primary_accession_ids, self._virus_db_id, self._sources)
                else:
                    current_sequences_in_db = database_tom.try_py_function(
                        vcm.sequence_alternative_accession_ids, self._virus_db_id, self._sources)
                current_sequences_in_db = Counter(current_sequences_in_db)
                inserted_in_db = current_sequences_in_db - self._sequences_in_db_at_start
                removed_from_db = self._sequences_in_db_at_start - current_sequences_in_db
                message += f"Sequences in DB:\n" \
                           f"\tprevious number {sum(self._sequences_in_db_at_start.values())}.\n" \
                           f"\tcurrent number {sum(current_sequences_in_db.values())}\n" \
                           f"\tdifference: {sum(inserted_in_db.values())} new - {sum(removed_from_db.values())} missing from source and deleted  = {sum(inserted_in_db.values()) - sum(removed_from_db.values())}\n"
                # check duplicates
                num_duplicates_at_start = sum(self._sequences_in_db_at_start.values()) - len(set(self._sequences_in_db_at_start))
                num_current_duplicates = sum(current_sequences_in_db.values()) - len(set(current_sequences_in_db))
                if num_duplicates_at_start > 0 or num_current_duplicates > 0:
                    warn = True
                    message += f"Duplicated accession_ids in DB:\n" \
                               f"\tprevious number of duplicates: {num_duplicates_at_start}\n" \
                               f"\tcurrent number of duplicates: {num_current_duplicates}\n" \
                               f"\tdetail of current duplicated accession_ids: {sorted(list(current_sequences_in_db - Counter(set(current_sequences_in_db))))}\n"

                # check errors in insertions
                changed_in_db = Counter(removed_samples_acc_id) # this set cannot be retrieved from the database
                completed_not_inserted = Counter(completed_samples_acc_id) - inserted_in_db - changed_in_db
                inserted_not_completed = inserted_in_db - Counter(completed_samples_acc_id) - changed_in_db
                if sum(completed_not_inserted.values()) > 0:
                    warn = True
                    message += f'Number of samples completed and not imported into the DB: {sum(completed_not_inserted.values())}\n' \
                               f'\taccession ids: {sorted(list(completed_not_inserted.elements()))}\n'
                if sum(inserted_not_completed.values()) < 0:
                    warn = True
                    message += f'Number of samples imported into the DB but incomplete: {sum(inserted_not_completed.values())}\n' \
                               f'\taccession ids: {sorted(list(inserted_not_completed.elements()))}\n'

            # performance
            message += self._performance_message(len(completed_samples_acc_id))

            if not warn:
                logger.info(message)
            else:
                logger.warning(message)

    @staticmethod
    def _queue_to_set(q: Queue):
        s = set()
        while not q.empty():
            try:
                s.add(q.get(False))
            except Empty:
                continue
        return s


_stats_type: Optional[StatsType] = None


def schedule_samples(set_info: StatsType):
    global _stats_type
    _stats_type = set_info


def completed_sample(*args, **kwargs):
    _stats_type.completed_sample(*args, **kwargs)


def removed_samples(*args, **kwargs):
    _stats_type.removed_samples(*args, **kwargs)


def check_samples_imported(*args, **kwargs):
    if _stats_type is not None:
        _stats_type.check_samples_imported(*args, **kwargs)
    else:
        logger.warning('\n'
                       'STATS MODULE:\n'
                       'Statistics not available because the program ended before this module was initialized.')
