import os
from enum import Enum

download_folder = f".{os.path.sep}downloads{os.path.sep}"
local_folder = download_folder + f"entrez_api{os.path.sep}"
local_folder_taxonomy = download_folder + f"entrez_api_taxonomy{os.path.sep}"
local_folder_nuc_variant_and_effects = download_folder + f"nuc_variants_and_effects{os.path.sep}"
local_folder_annotations_and_aa_var = download_folder + f"annotations_and_aa_variants{os.path.sep}"
local_folder_sequence = download_folder + f"sequence{os.path.sep}"

base_folder = f".{os.path.sep}generated{os.path.sep}"

class FileType(Enum):
    TaxonomyData = 'taxonomy'
    SequenceOrSampleData = 'samples'
    NucleotideVariants = 'nuc_variants'
    Logs = 'logs'

def get_local_folder_for(source_name: str, _type: FileType) -> str:
    """
    :param source_name: the name of the source requiring the directory. As this name will be part of the path, it is
    important that a source issues always the same name in order to get a stable path.
    :param _type: one of locations.FileType
    :return: the path to an existing and usable folder.
    """
    path = f"{base_folder}{source_name}{os.path.sep}{_type.value}{os.path.sep}"
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    return path

