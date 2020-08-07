from code.utils import create_dir
from code.epitopes.BepiPredEpitopes import BepiPredEpitopes

###               SARS-CoV-2                 ###
### Bepipred 2.0 output  --> epitopes csv    ###
###                                          ###

#input
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/sars_cov2_data/pred_epitopes_input"
pred_epitopes_folder = "bepipred_out"
pred_file = "bepipred_sars_cov2.tsv"
virus_ncbi_id = "NC_045512.2"

# output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/pred_epitopes/bepipred"
create_dir(out_path)

BepiPredEpitopes_sars_cov2 = BepiPredEpitopes(data_path,pred_epitopes_folder,virus_ncbi_id,out_path)

BepiPredEpitopes_sars_cov2.save2tsv(pred_file)