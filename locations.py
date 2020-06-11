import os

download_folder = "downloads/"
local_folder = download_folder + "entrez_api/"
local_folder_taxonomy = download_folder + "entrez_api_taxonomy/"
local_folder_variant = download_folder + "variant/"
local_folder_sequence = download_folder +"sequence/"


def create_local_folders():
    for folder in [local_folder, local_folder_taxonomy, local_folder_variant, local_folder_sequence]:
        if not os.path.exists(folder):
            os.makedirs(folder)
