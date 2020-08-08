from VirusGenoUtil.code.utils import create_dir
from VirusGenoUtil.code.epitopes.IEDBEpitopes import IEDBEpitopes

###                                                          ###
###                  Import IEDB epitopes                    ###
###                                                          ###

###                                                          ###
###     input = IEDB Bcells.csv + Tcells.csv                 ###
###     output = epitope csv for all taxon ids               ###
###                                                          ###

# input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/iedb_data"
cells_epitopes_folder = "cell_epitopes"
virus_proteins_folder = "viruses_proteins"
host_taxon, host_name = "taxon_9606", "Homo Sapiens"
assay_type = "positive"

# output
output_path = "/home/damian/Documents/L3S/projects/sars_cov2/IEDB_epitopes"
create_dir(output_path)

IEDBEpitopes = IEDBEpitopes(data_path, cells_epitopes_folder, virus_proteins_folder,
                            host_taxon, host_name, assay_type, output_path)
IEDBEpitopes.process_all_viruses()
