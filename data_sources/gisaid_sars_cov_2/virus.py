import json
import os

from json import JSONDecodeError
from typing import Optional

from loguru import logger
from tqdm import tqdm

from data_sources.gisaid_sars_cov_2.sample import GISAIDSarsCov2Sample
from data_sources.virus import VirusSource
import stats_module


# noinspection PyMethodMayBeStatic
from locations import get_local_folder_for, FileType


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

    def virus_samples(self, virus_id: int, from_sample: Optional[int] = None, to_sample: Optional[int] = None):
        with open(self.data_path, mode='r') as input_file:
            num_lines = sum(1 for line in input_file)
            input_file.seek(0, 0)   # reset pointer
            lines_read = 0
            lines_to_read = num_lines if (from_sample is None or to_sample is None) else to_sample - from_sample
            progress = tqdm(total=lines_to_read)
            stats_module.schedule_samples(stats_module.StatsBasedOnTotals(lines_to_read, virus_id, ['GISAID']))
            for line in input_file:
                if from_sample is not None and to_sample is not None:
                    if lines_read < from_sample:
                        lines_read += 1
                        continue    # read line without action
                    elif lines_read >= to_sample:
                        return      # terminate loop
                try:
                    lines_read += 1
                    progress.update()
                    yield GISAIDSarsCov2Sample(json.loads(line))
                except JSONDecodeError:
                    pass

