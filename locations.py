import os

download_folder = f".{os.path.sep}/downloads{os.path.sep}"
local_folder = download_folder + f"entrez_api{os.path.sep}"
local_folder_taxonomy = download_folder + f"entrez_api_taxonomy{os.path.sep}"
local_folder_nuc_variant_and_effects = download_folder + f"nuc_variants_and_effects{os.path.sep}"
local_folder_annotations_and_aa_var = download_folder + f"annotations_and_aa_variants{os.path.sep}"
local_folder_sequence = download_folder + f"sequence{os.path.sep}"


def create_local_folders():
    for folder in [local_folder, local_folder_taxonomy, local_folder_nuc_variant_and_effects, local_folder_sequence]:
        if not os.path.exists(folder):
            os.makedirs(folder)
