"""
Created by tomalf2 on nov, 2020.
"""
from jsonschema import validate, ValidationError, SchemaError
import json
from os.path import sep
from loguru import logger


_db_config_schema_path = f".{sep}db_config{sep}db_import_configuration.schema.json"
_db_config_path = f".{sep}db_config{sep}db_import_configuration.json"
_args_are_valid: bool = False

_db_name = None # this is received from program arguments and set with "set_db_name()"


def _validate_import_configuration_file():
    with open(_db_config_path) as json_file:
        configuration: dict = json.load(json_file)
    with open(_db_config_schema_path) as json_schema_file:
        configuration_schema: dict = json.load(json_schema_file)
    try:
        validate(configuration, configuration_schema)
    except SchemaError as e:
        logger.error("The JSON schema file is syntactically valid but contains a semantic error.")
        raise e
    except ValidationError as e:
        details = f"""DETAILS OF ERROR:\n---------\n{e.absolute_path}\n---------\n{e.absolute_schema_path}"""
        logger.error("the provided configuration file does not validate against the schema. " + details)
        raise e
    else:
        return True


def get_database_config_params():
    global _args_are_valid
    if not _args_are_valid:
        _args_are_valid = _validate_import_configuration_file()

    with open(_db_config_path) as json_file:
        configuration: dict = json.load(json_file)
    configuration["db_name"] = _db_name
    return configuration


def set_db_name(db_name):
    global _db_name
    _db_name = db_name