"""
Created by tomalf2 on nov, 2020.
"""
from loguru import logger
import json
from os.path import sep
import time
import pprint


_overlaps_config_path = f".{sep}db_config{sep}db_overlaps_configuration.json"
_read_first_time = True

def get_import_params_for(key: str) -> dict:
    with open(_overlaps_config_path) as json_file:
        configuration: dict = json.load(json_file)
    global _read_first_time
    if _read_first_time:
        logger.warning("DBs configuration for overlap checks:\n" + pprint.pformat(configuration) +
                       "\nCheck database names are correct. The program will resume in 15s.")
        _read_first_time = False
        time.sleep(15)
    try:
        source_params = configuration[key]
    except KeyError:
        known_importable_sources = configuration.keys()
        sources_list =  '\n-'.join(known_importable_sources)
        raise KeyError(f"No import configuration found for source named {key}\n"
                       f"Available importable sources are:"
                       f"{sources_list}")
    return source_params
