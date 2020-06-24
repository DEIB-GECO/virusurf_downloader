import os

download_folder = f".{os.path.sep}/downloads{os.path.sep}"
local_folder = download_folder + f"entrez_api{os.path.sep}"
local_folder_taxonomy = download_folder + f"entrez_api_taxonomy{os.path.sep}"
local_folder_variant = download_folder + f"variant{os.path.sep}"
local_folder_sequence = download_folder + f"sequence{os.path.sep}"


def create_local_folders():
    for folder in [local_folder, local_folder_taxonomy, local_folder_variant, local_folder_sequence]:
        if not os.path.exists(folder):
            os.makedirs(folder)
