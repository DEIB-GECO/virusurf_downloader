from code.utils import create_dir
from code.epitopes.IEDBEpitopes import IEDBEpitopes
from pathlib import Path
from os.path import join
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
root = Path.cwd().parent
data_path = join(root,"data/iedb_input")
cells_epitopes_folder = "cell_epitopes"
virus_proteins_folder = "viruses_proteins"
host_taxon, host_name = "taxon_9606", "Homo Sapiens"
assay_type = "positive"

# output
output_path = join(root,"data/iedb_output")

IEDBEpitopes = IEDBEpitopes(data_path, cells_epitopes_folder, virus_proteins_folder,
                            host_taxon, host_name, assay_type, output_path)

# extract IEDB epitopes from virus:
virus_taxid = 2697049
IEDBEpitopes.process_virus(virus_taxid)

# then get epitopes and epitope fragments
virus_epitopes = IEDBEpitopes.get_virus_epitopes()
print(virus_epitopes[0])
virus_epi_fragments = IEDBEpitopes.get_virus_epi_fragments()
print(virus_epi_fragments[0])