import os, sys
sys.path.insert(0, os.path.abspath(".."))

from VirusGenoUtil.code.utils import create_dir
from VirusGenoUtil.code.epitopes.IEDBEpitopes import IEDBEpitopes
from pathlib import Path
from os.path import join
from os.path import sep
###                                                          ###
###           Import IEDB epitopes in ViruSurf               ###
###                                                          ###

###                                                          ###
###     input = IEDB Bcells.csv + Tcells.csv                 ###
###     output = epitope csv for all taxon ids               ###
###                                                          ###

###             Paths                 ###
###  absolute paths for damianos pc   ###
###  relative paths for ViruSurf code ###

# input
data_path = f"./VirusGenoUtil{sep}data{sep}iedb_input"
cells_epitopes_folder = "cell_epitopes"
virus_proteins_folder = "viruses_proteins"
ontie_download_folder = "ONTIE_downloads"
host_taxon, host_name = "taxa_all", None
assay_type = "all"

# output
output_path = f"./VirusGenoUtil{sep}data{sep}iedb_output"


def epitopes_for_virus_taxon(vir_taxon_id: int):
    # extract IEDB epitopes from virus:
    epitopes = IEDBEpitopes(data_path, cells_epitopes_folder, virus_proteins_folder, ontie_download_folder,
                            host_taxon, host_name, assay_type, output_path)
    epitopes.process_virus(vir_taxon_id)
    # get epitopes and epitope fragments
    return epitopes.get_virus_epitopes(), epitopes.get_virus_epi_fragments()
