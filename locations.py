import os
from enum import Enum

from loguru import logger

_generated_folder = f".{os.path.sep}generated{os.path.sep}"


class FileType(Enum):
    TaxonomyData = 'taxonomy'
    SequenceOrSampleData = 'samples'
    NucleotideVariants = 'nuc_variants'
    Logs = 'logs'
    Annotations = 'annotations'
    Fasta = 'fasta'
    HostData = 'hosts'


def get_local_folder_for(source_name: str, _type: FileType) -> str:
    """
    :param source_name: the name of the source requiring the directory. As this name will be part of the path, it is
    important that a source issues always the same name in order to get a stable path.
    :param _type: one of locations.FileType
    :return: the path to an existing and usable folder.
    """
    path = f"{_generated_folder}{source_name}{os.path.sep}{_type.value}{os.path.sep}"
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    return path


def remove_file(file_path):
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except OSError as e:
            logger.error(f"Failed to remove file {file_path} with error: {e.strerror}")
