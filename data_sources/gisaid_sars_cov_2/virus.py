import json
import os

from json import JSONDecodeError
from typing import Generator

from loguru import logger
from tqdm import tqdm

from data_sources.gisaid_sars_cov_2.sample import GISAIDSarsCov2Sample
from data_sources.virus import VirusSource
from db_config.database_tom import try_py_function, Sequence, HostSample, SequencingProject


# noinspection PyMethodMayBeStatic


class GISAIDSarsCov2(VirusSource):

    name = 'GISAID_sars_cov_2'
    data_path = f'.{os.path.sep}..{os.path.sep}GISAID{os.path.sep}export.json'

    def __init__(self):
        super().__init__()
        logger.info(f'importing virus {GISAIDSarsCov2.name}')

    def taxon_id(self):
        return 2697049

    def taxon_name(self):
        return 'Severe acute respiratory syndrome coronavirus 2'

    def family(self):
        return 'Coronaviridae'

    def sub_family(self):
        return 'Orthocoronavirinae'

    def genus(self):
        return 'Betacoronavirus'

    def species(self):
        return 'Severe acute respiratory syndrome-related coronavirus'

    def equivalent_names(self):
        return 'SARS-CoV-2, 2019-nCoV, COVID-19, COVID-19 virus, COVID19, HCoV-19, Human coronavirus 2019, SARS-2, SARS-CoV2, SARS2, Wuhan coronavirus, Wuhan seafood market pneumonia virus'

    def molecule_type(self):
        return 'RNA'

    def is_single_stranded(self):
        return True

    def is_positive_stranded(self):
        return True

    def count_sequences_in_file(self):
        with open(self.data_path, mode='r') as input_file:
            num_lines = sum(1 for line in input_file)
        return num_lines

    # def virus_samples(self, virus_id: int, from_sample: Optional[int] = None, to_sample: Optional[int] = None):
    #     num_lines = self.count_sequences_in_file()
    #     with open(self.data_path, mode='r') as input_file:
    #         lines_read = 0
    #         lines_to_read = num_lines if (from_sample is None or to_sample is None) else to_sample - from_sample
    #         progress = tqdm(total=lines_to_read)
    #         stats_module.schedule_samples(stats_module.StatsBasedOnTotals(lines_to_read, virus_id, ['GISAID']))
    #         for line in input_file:
    #             if from_sample is not None and to_sample is not None:
    #                 if lines_read < from_sample:
    #                     lines_read += 1
    #                     continue    # read line without action
    #                 elif lines_read >= to_sample:
    #                     return      # terminate loop
    #             try:
    #                 lines_read += 1
    #                 progress.update()
    #                 yield GISAIDSarsCov2Sample(json.loads(line))
    #             except JSONDecodeError:
    #                 pass

    def get_sequences_in_current_data(self) -> dict:
        def do(session):
            all_sequences = dict()
            for db_items in session.query(Sequence, HostSample, SequencingProject) \
                                      .filter(Sequence.strain_name.isnot(None),
                                              Sequence.host_sample_id == HostSample.host_sample_id,
                                              Sequence.sequencing_project_id == SequencingProject.sequencing_project_id) \
                                      .yield_per(100):
                source_seq, source_host, source_prj = db_items
                all_sequences[source_seq.accession_id] = (source_seq.strain_name, source_seq.length,
                                                          source_seq.gc_percentage, source_seq.n_percentage,
                                                          source_host.collection_date, source_host.originating_lab,
                                                          source_prj.submission_date, source_prj.sequencing_lab,
                                                          source_host.country, source_host.region,
                                                          source_host.isolation_source)
            return all_sequences

        return try_py_function(do)

    def get_sequences_of_updated_source(self) -> Generator[GISAIDSarsCov2Sample, None, None]:
        with open(self.data_path, mode='r') as input_file:
            for line in input_file:
                try:
                    yield GISAIDSarsCov2Sample(json.loads(line))
                except JSONDecodeError:
                    pass

    def deltas(self):
        """
        Returns the set of Sequence.accession_id to remove from current database, and the set of Sequence.accession_id
        that should be imported from the source. Note that because sequences in the current data may have been updated,
        the two sets are not necessarily disjoint. So first remove the sequences to be removed and only then
         insert the new ones.
         """
        acc_id_remote = set([x.primary_accession_number() for x in self.get_sequences_of_updated_source()])  # all the sequences from remote
        logger.info('Collecting current sequences...')
        current_data = self.get_sequences_in_current_data()
        acc_id_current = set(current_data.keys())

        acc_id_missing_in_remote = acc_id_current - acc_id_remote
        acc_id_missing_in_current = acc_id_remote - acc_id_current

        acc_id_present_in_current_and_remote = acc_id_current & acc_id_remote
        acc_id_changed = []

        # to compare sequences present in both, have to scan the file of remote sequences
        logger.info('Comparing current data with updated source data')
        for new_sequence in tqdm(self.get_sequences_of_updated_source()):
            acc_id = new_sequence.primary_accession_number()
            try:
                current_sequence_data = current_data[acc_id]
                # if the gisaid id is the same, compare metadata to detect if the sample changed over time
                if current_sequence_data[0] != new_sequence.strain() \
                        or current_sequence_data[1] != new_sequence.length() \
                        or current_sequence_data[2] != new_sequence.gc_percent() \
                        or current_sequence_data[3] != new_sequence.n_percent() \
                        or current_sequence_data[4] != str(new_sequence.collection_date()) \
                        or current_sequence_data[5] != new_sequence.originating_lab() \
                        or str(current_sequence_data[6]) != str(new_sequence.submission_date()) \
                        or current_sequence_data[7] != new_sequence.sequencing_lab() \
                        or (current_sequence_data[8], current_sequence_data[9]) != new_sequence.country__region__geo_group()[:2] \
                        or current_sequence_data[9] != new_sequence.isolation_source():
                    acc_id_changed.append(acc_id)
            except KeyError:
                pass  # the accession id is not present in current data. it's a new sequence

        # compute additional sets
        acc_id_changed = set(acc_id_changed)
        acc_id_unchanged = acc_id_present_in_current_and_remote - acc_id_changed

        acc_id_to_remove = acc_id_missing_in_remote | acc_id_changed
        acc_id_to_import = acc_id_missing_in_current | acc_id_changed

        logger.info(f'\n'
                    f'# Sequences from remote source: {len(acc_id_remote)}. Of which\n'
                    f'# {len(acc_id_unchanged)} present also locally and unchanged\n'
                    f'# {len(acc_id_missing_in_current)} never seen before.\n'
                    f'# {len(acc_id_changed)} have different attributes in the remote source\n'
                    f'# Sequences from local source: {len(acc_id_current)}. Of which\n'
                    f'# {len(acc_id_missing_in_remote)} are missing from remote and must be removed from local.\n'
                    f'# In conclusion: {len(acc_id_to_remove)} sequences will be removed because missing or changed in remote\n'
                    f'# and {len(acc_id_to_import)} sequences will be imported because novel or changed in remote.')

        return acc_id_to_remove, acc_id_to_import

