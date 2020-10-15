from VirusGenoUtil.code.variants.MFA2CSV import MFA2CSV

### ### ### ### ### ### ### ### ### ### ### ####
# virus reference + target genomes -> variants #
### ### ### ### ### ### ### ### ### ### ### ####

###
# install virulign
###

###
# run the provided bash for your specified (target) sequences
###

###
# convert MFA --> variants tabular file per target
###

###                    ###
###     SARS-Cov2      ###
###                    ###

# ===========#
#   Test     #
# ===========#
# use as test orf the ORF7b for SARS-cov2 from virulign and create 3 sequences with the following variants:
# >Ref
# CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAA
# >Seq1
# ------------ATTGACACTGTACCAGTAACATTAAAGCCAGGAATGGATGGACCAAAG
# >Seq2
# CCTATG---------GAAACTGTGCCAGTAAAATTAAAGCCAGGAATGGAT--------- (not displayed)
# >Seq3
# CTCATTAGTCCTATTAGT---------GTAAAATTAAAACCAGGAATGGATGGCCCAAGG  (not displayed)
# >Seq4
# ------AGTCCCATTGAAACTGTACCAGTAAAA---------GGA---GATGGCCCAAAG (not displayed)
"""
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
alignments_folder = "test_alignments"
xmls_folder = "test_orf_xml"
#ORF input:
mfa_file = "test_mfa_nt.fasta"
orf_xml_file = "test_orf.xml"
ncbi_ref_id = "test_ref"
#output
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/test_variants"
"""

"""
#======================================================#
# Example for S (Spike) protein and 2 target sequences #
#======================================================#
# inputs:
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
alignments_folder = "alignments"
xmls_folder = "xmls"
# ORF input:
mfa_file = "S_mfa_nt.fasta"
orf_xml_file = "S.xml"
ncbi_ref_id = "NC_045512.2"

# output:
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/variants"

print("==== ====")
MFA2CSV_S_protein_test = MFA2CSV(data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path)
MFA2CSV_S_protein_test.run(orf_xml_file, mfa_file)
print("=== * ===")
print("== *** ==")
"""

"""
# =================================== #
# Convert MFA to variants csv file    #
# for all ORFS of 10 target sequences #
# =================================== #

# input #
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
alignments_folder = "alignments"
xmls_folder = "xmls"
ncbi_ref_id = "NC_045512.2"
# output #
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/variants"
print("==== ====")
MFA2CSV_all_orf_10targets = MFA2CSV(data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path)
MFA2CSV_all_orf_10targets.run_multiple_orfs()
print("=== * ===")
print("== *** ==")
"""

# =================================== #
# Convert MFA to variants csv file    #
# for all ORFS of NCBI sequences      #
# =================================== #
# input #
data_path = "/home/damian/Documents/L3S/projects/sars_cov2/data"
alignments_folder = "ncbi_alignments"
xmls_folder = "xmls"
ncbi_ref_id = "NC_045512.2"
is_amino_acid = True
# output #
out_path = "/home/damian/Documents/L3S/projects/sars_cov2/ncbi_variants"
print("==== ====")
MFA2CSV_all_orf_ncbi_target_aa = MFA2CSV(data_path, alignments_folder, xmls_folder, ncbi_ref_id, out_path)
MFA2CSV_all_orf_ncbi_target_aa.run_multiple_orfs(is_amino_acid)
print("=== * ===")
print("== *** ==")
