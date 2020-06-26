import json

from json import JSONDecodeError
from loguru import logger
from tqdm import tqdm

from data_sources.gisaid_sars_cov_2.sample import GISAIDSarsCov2Sample
from data_sources.virus import VirusSource


# noinspection PyMethodMayBeStatic
class GISAIDSarsCov2(VirusSource):

    internal_name = 'GISAID_sars_cov_2'
    data_path = './downloads/gisaid_test/export.json'

    def __init__(self):
        super().__init__()
        logger.info(f'importing virus {GISAIDSarsCov2.internal_name}')

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

    def virus_samples(self):
        with open(GISAIDSarsCov2.data_path, mode='r') as input_file:
            for line in tqdm(input_file):
                try:
                    yield GISAIDSarsCov2Sample(json.loads(line))
                except JSONDecodeError:
                    pass

