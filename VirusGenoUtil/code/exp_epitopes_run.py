from code.epitopes.Immunodominance import Immunodominance
from code.epitopes.HomologyBasedEpitopes import HomologyBasedEpitopes

from code.utils import create_dir

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
# 1)sars cov1 proteins + immunome browser -> lower bound sliding window per position#
# 2)sars cov2 proteins -BLAST-> sars cov1 -> similarity per position                #
# 3)lower bound + similarity per position -> experimental epitopes for sars cov2    #
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #


###                                          ###
### B cells epitopes sars-cov1 --> sars-cov2 ###
###                                          ###

# input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/sars_cov2_data/exp_epitopes_input"
query_proteins_folder = "sars_cov1_proteins_Bcells"
exp_epitopes_folder = "Bcells"
assay_species = "human"
exp_epitopes_folder = exp_epitopes_folder
window_size = 10
save_plot = True
sliding_avg_cutoff = 0.3
split_at_min_RF = True
save_immunodominant_reg = True
process_tcells = False

# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Bcells"
create_dir(out_path)


"""
print("==== ====")
Immunodominance_Bcells = Immunodominance(data_path, query_proteins_folder, exp_epitopes_folder, assay_species, out_path)
processed_proteins = Immunodominance_Bcells.process_all_proteins(process_tcells, window_size, save_plot, sliding_avg_cutoff,
                                                                 split_at_min_RF,
                                                                 save_immunodominant_reg)
"""


print("==== ====")
# input
subject_proteins_folder = "sars_cov2_proteins"
blast_folder = "blast"
ncbi_ids = {"relative_organism": "NC_004718.3", "target_organism": "NC_045512.2"}

HomologyBasedEpitopes_Bcells = HomologyBasedEpitopes(ncbi_ids, data_path, subject_proteins_folder, blast_folder,
                                                     out_path)
HomologyBasedEpitopes_Bcells.find_epitopes(process_tcells)


###                                          ###
### T cells epitopes sars-cov1 --> sars-cov2 ###
###                                          ###
"""
# input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/sars_cov2_data/exp_epitopes_input"
query_proteins_folder = "sars_cov1_proteins_Tcells"
exp_epitopes_folder = "Tcells"
assay_species = "human"
window_size = 10
save_plot = True
sliding_avg_cutoff = 0.3
split_at_min_RF = True
save_immunodominant_reg = True
process_tcells = True

# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/exp_epitopes/Tcells"
create_dir(out_path)
"""

"""
print("==== ====")
Immunodominance_Tcells = Immunodominance(data_path, query_proteins_folder, exp_epitopes_folder, assay_species, out_path)
processed_proteins = Immunodominance_Tcells.process_all_proteins(process_tcells, window_size, save_plot,
                                                                 sliding_avg_cutoff, split_at_min_RF,
                                                                 save_immunodominant_reg)
"""

"""
print("==== ====")
# input
subject_proteins_folder = "sars_cov2_proteins"
blast_folder = "blast"
ncbi_ids = {"relative_organism": "NC_004718.3", "target_organism": "NC_045512.2"}

HomologyBasedEpitopes_Tcells = HomologyBasedEpitopes(ncbi_ids, data_path, subject_proteins_folder, blast_folder,
                                                     out_path)
HomologyBasedEpitopes_Tcells.find_epitopes(process_tcells)
"""

print("=== * ===")
print("== *** ==")
