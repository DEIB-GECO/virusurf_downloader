#!/bin/bash

# make sure this script is executable:  chmod 755 myscript.sh
# set your user as owner of the script: chown myusername: myscript.sh

echo "Assumptions for the proper execution of this script:
- command line arguments: <virusurf_directory_path> <database_name> <download_updated_epitopes_?>;
- the postgres password file contains one entry for the given database formatted as host:port:dat_name:user:psw;
- a conda environment named 'vcm' exists and it has the packages necessary for virusurf_downlaoder already installed;"
sleep 10

# input parameters
virusurf_dir=${1}/
database_name=$2
enable_download=$(echo "$3" | tr '[:upper:]' '[:lower:]') # convert lowercase
echo "Check of arguments:
- VIRUSURF DIR: ${virusurf_dir};
- VIRUSURF DB: ${database_name};
- DOWNLOAD: $(if [ "$enable_download" == "yes" ] || [ "$enable_download" == "true" ]; then echo true; else echo false; fi);
The program resumes in 10 seconds."
sleep 10

log_file_path=${virusurf_dir}logs/epitope_others_updates.log
notifier_script_path=${virusurf_dir}bash_scripts/notify_me.sh

# define utility functions
timestamp() {
  TZ=Europe/Rome date +"%Y_%m_%d__%H_%M_%S" # current time
}
notify_me() { # expects a message as argument
  if test -f "$notifier_script_path"; then
    /bin/bash "$notifier_script_path" "$1"
  else
    echo "virusurf self-update notification not available"  | tee -a $log_file_path
  fi
}
check_exit_code() {    # expects exit status code as argument
  if [ "$1" -ne 0 ]; then
    echo "* Last command terminated abnormally (exit code ${1}). Script interrupted." | tee -a $log_file_path
    notify_me "* $(timestamp) virusurf self-update notification: one of the update steps failed. Check $log_file_path
    for details." | tee -a $log_file_path
    echo ""
    echo ""
    echo ""
    exit 1
  fi
}

# switch to conda environment vcm
echo "Switching to  conda environment 'vcm'."
source ~/anaconda3_new/etc/profile.d/conda.sh  # finds command conda
check_exit_code "$?"
conda activate vcm
check_exit_code "$?"

##############################################################################

echo "The output of this script is appended to ${log_file_path}"

# Begin logging of operations
echo "################################          ##################################" | tee -a $log_file_path
echo "##############   Epitope ALL-but-SC2 Updater for ViruSurf   ################" | tee -a $log_file_path
echo "#######            - downloads updates from IEDB                     #######" | tee -a $log_file_path
echo "#######            - deletes outdated ALL-but-SC2 epitopes           #######" | tee -a $log_file_path
echo "#######            - import ALL-but-SC2-related epitopes             #######" | tee -a $log_file_path
echo "#######            - refresh ALL-but-SC2-related views for EpiSurf   #######" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic epitopes update script started at $(timestamp)" | tee -a $log_file_path


# Update epitope files from IEDB
if [ "$enable_download" == "yes" ] || [ "$enable_download" == "true" ]; then
  echo "* Download of updated source files from IEDB at $(timestamp) ..." | tee -a $log_file_path
  cd $virusurf_dir
  python main.py download epitopes
  check_exit_code "$?"
fi


# Delete old epitopes
echo "* Deletion of old non-SC2 epitopes at $(timestamp) ..." | tee -a $log_file_path
cd $virusurf_dir
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -c "
-- delete epitope fragments of this virus
WITH which_virus AS (
    SELECT virus_id
    FROM virus
    WHERE taxon_id != 2697049 -- exclude SC2 taxon_id
)
DELETE FROM epitope_fragment AS ef
USING epitope AS e
WHERE ef.epitope_id = e.epitope_id AND e.virus_id IN (SELECT which_virus.virus_id FROM which_virus);

-- delete epitopes of this virus
WITH which_virus AS (
    SELECT virus_id
    FROM virus
    WHERE taxon_id != 2697049 -- exclude SC2 taxon_id
)
DELETE FROM epitope AS e
WHERE e.virus_id IN (SELECT which_virus.virus_id FROM which_virus);
" | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"

# Insert ALL-but-SC2 epitopes
cd $virusurf_dir
echo "* Begin import of epitopes for virus bundibugyo_ebolavirus at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} bundibugyo_ebolavirus
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_1 at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_1
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_2 at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_2
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_3 at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_3
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_4 at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_4
check_exit_code "$?"

echo "* Begin import of epitopes for virus mers at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} mers
check_exit_code "$?"

echo "* Begin import of epitopes for virus sars_cov_1 at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} sars_cov_1
check_exit_code "$?"

echo "* Begin import of epitopes for virus sudan_ebolavirus at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} sudan_ebolavirus
check_exit_code "$?"

echo "* Begin import of epitopes for virus zaire_ebolavirus at $(timestamp) " | tee -a $log_file_path
python main.py epitopes ${database_name} zaire_ebolavirus
check_exit_code "$?"


# Refresh materialized views for episurf
echo "* Refresh ALL-but-SC2 eptiope variants and info materialized view at $(timestamp) "
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -c "
-- bombali ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_2010960_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_2010960_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_2010960_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_2010960_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_2010960_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_2010960_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_2010960_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_2010960_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_2010960_rna_depende;

-- mers
REFRESH MATERIALIZED VIEW public.epitope_1335626_1ab_polypro;

REFRESH MATERIALIZED VIEW public.epitope_1335626_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_1335626_hel;

REFRESH MATERIALIZED VIEW public.epitope_1335626_exon;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nendou;

REFRESH MATERIALIZED VIEW public.epitope_1335626_2_o_methylt;

REFRESH MATERIALIZED VIEW public.epitope_1335626_1a_polyprot;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp1_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp2_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp3_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp4_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp5_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp6_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp7_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp8_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp9_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp10_prote;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nsp11_prote;

REFRESH MATERIALIZED VIEW public.epitope_1335626_spike_prote;

REFRESH MATERIALIZED VIEW public.epitope_1335626_ns3_protein;

REFRESH MATERIALIZED VIEW public.epitope_1335626_ns4a_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_ns4b_protei;

REFRESH MATERIALIZED VIEW public.epitope_1335626_ns5_protein;

REFRESH MATERIALIZED VIEW public.epitope_1335626_envelope_pr;

REFRESH MATERIALIZED VIEW public.epitope_1335626_membrane_pr;

REFRESH MATERIALIZED VIEW public.epitope_1335626_nucleocapsi;

REFRESH MATERIALIZED VIEW public.epitope_1335626_orf8b_prote;

-- sars cov 1
REFRESH MATERIALIZED VIEW public.epitope_694009_orf1ab_poly;

REFRESH MATERIALIZED VIEW public.epitope_694009_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_694009_helicase_nt;

REFRESH MATERIALIZED VIEW public.epitope_694009_3_to_5_exon;

REFRESH MATERIALIZED VIEW public.epitope_694009_endoribonuc;

REFRESH MATERIALIZED VIEW public.epitope_694009_2_o_mtase;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf1a_polyp;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp1;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp2;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp3;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp4;

REFRESH MATERIALIZED VIEW public.epitope_694009_3c_like_pro;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp6;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp7;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp8;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp9;

REFRESH MATERIALIZED VIEW public.epitope_694009_nsp10;

REFRESH MATERIALIZED VIEW public.epitope_694009_ndp11;

REFRESH MATERIALIZED VIEW public.epitope_694009_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf3a_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf3b_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_small_envel;

REFRESH MATERIALIZED VIEW public.epitope_694009_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf6_protei;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf7a_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf7b_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf8a_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf8b_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_nucleocapsi;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf9b_prote;

REFRESH MATERIALIZED VIEW public.epitope_694009_orf9a_prote;

-- tai fores ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_186541_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186541_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186541_vp35;

REFRESH MATERIALIZED VIEW public.epitope_186541_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_186541_vp40;

REFRESH MATERIALIZED VIEW public.epitope_186541_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_186541_sgp;

REFRESH MATERIALIZED VIEW public.epitope_186541_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_186541_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_186541_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_186541_vp30;

REFRESH MATERIALIZED VIEW public.epitope_186541_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_186541_vp24;

REFRESH MATERIALIZED VIEW public.epitope_186541_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_186541_polymerase;

REFRESH MATERIALIZED VIEW public.epitope_186541_rna_depende;

-- bundibugyo ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_565995_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_565995_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_565995_vp35;

REFRESH MATERIALIZED VIEW public.epitope_565995_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_565995_vp40;

REFRESH MATERIALIZED VIEW public.epitope_565995_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_565995_sgp;

REFRESH MATERIALIZED VIEW public.epitope_565995_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_565995_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_565995_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_565995_vp30;

REFRESH MATERIALIZED VIEW public.epitope_565995_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_565995_vp24;

REFRESH MATERIALIZED VIEW public.epitope_565995_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_565995_polymerase;

REFRESH MATERIALIZED VIEW public.epitope_565995_rna_depende;

-- reston ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_186539_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186539_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186539_vp35;

REFRESH MATERIALIZED VIEW public.epitope_186539_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_186539_vp40;

REFRESH MATERIALIZED VIEW public.epitope_186539_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_186539_sgp;

REFRESH MATERIALIZED VIEW public.epitope_186539_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_186539_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_186539_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_186539_vp30;

REFRESH MATERIALIZED VIEW public.epitope_186539_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_186539_vp24;

REFRESH MATERIALIZED VIEW public.epitope_186539_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_186539_polymerase;

REFRESH MATERIALIZED VIEW public.epitope_186539_rna_depende;

-- sudan ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_186540_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186540_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186540_vp35;

REFRESH MATERIALIZED VIEW public.epitope_186540_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_186540_vp40;

REFRESH MATERIALIZED VIEW public.epitope_186540_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_186540_sgp;

REFRESH MATERIALIZED VIEW public.epitope_186540_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_186540_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_186540_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_186540_vp30;

REFRESH MATERIALIZED VIEW public.epitope_186540_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_186540_vp24;

REFRESH MATERIALIZED VIEW public.epitope_186540_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_186540_polymerase;

REFRESH MATERIALIZED VIEW public.epitope_186540_rna_depende;

-- zaire ebolavirus
REFRESH MATERIALIZED VIEW public.epitope_186538_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186538_nucleoprote;

REFRESH MATERIALIZED VIEW public.epitope_186538_vp35;

REFRESH MATERIALIZED VIEW public.epitope_186538_polymerase_;

REFRESH MATERIALIZED VIEW public.epitope_186538_vp40;

REFRESH MATERIALIZED VIEW public.epitope_186538_matrix_prot;

REFRESH MATERIALIZED VIEW public.epitope_186538_sgp;

REFRESH MATERIALIZED VIEW public.epitope_186538_spike_glyco;

REFRESH MATERIALIZED VIEW public.epitope_186538_small_secre;

REFRESH MATERIALIZED VIEW public.epitope_186538_second_secr;

REFRESH MATERIALIZED VIEW public.epitope_186538_vp30;

REFRESH MATERIALIZED VIEW public.epitope_186538_minor_nucle;

REFRESH MATERIALIZED VIEW public.epitope_186538_vp24;

REFRESH MATERIALIZED VIEW public.epitope_186538_membrane_as;

REFRESH MATERIALIZED VIEW public.epitope_186538_polymerase;

REFRESH MATERIALIZED VIEW public.epitope_186538_rna_depende;

-- dengue virus 4
REFRESH MATERIALIZED VIEW public.epitope_11070_polyprotein;

REFRESH MATERIALIZED VIEW public.epitope_11070_anchored_ca;

REFRESH MATERIALIZED VIEW public.epitope_11070_capsid_prot;

REFRESH MATERIALIZED VIEW public.epitope_11070_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11070_protein_pr;

REFRESH MATERIALIZED VIEW public.epitope_11070_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11070_envelope_pr;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_protein_2k;

REFRESH MATERIALIZED VIEW public.epitope_11070_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11070_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_11070_sfrna2;

REFRESH MATERIALIZED VIEW public.epitope_11070_sfrna3;

REFRESH MATERIALIZED VIEW public.epitope_11070_sfrna4;

-- dengue virus 3
REFRESH MATERIALIZED VIEW public.epitope_11069_polyprotein;

REFRESH MATERIALIZED VIEW public.epitope_11069_anchored_ca;

REFRESH MATERIALIZED VIEW public.epitope_11069_capsid_prot;

REFRESH MATERIALIZED VIEW public.epitope_11069_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11069_protein_pr;

REFRESH MATERIALIZED VIEW public.epitope_11069_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11069_envelope_pr;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_protein_2k;

REFRESH MATERIALIZED VIEW public.epitope_11069_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11069_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_11069_sfrna1;

REFRESH MATERIALIZED VIEW public.epitope_11069_sfrna2;

REFRESH MATERIALIZED VIEW public.epitope_11069_sfrna3;

REFRESH MATERIALIZED VIEW public.epitope_11069_sfrna4;

-- dengue virus 2
REFRESH MATERIALIZED VIEW public.epitope_11060_polyprotein;

REFRESH MATERIALIZED VIEW public.epitope_11060_anchored_ca;

REFRESH MATERIALIZED VIEW public.epitope_11060_capsid_prot;

REFRESH MATERIALIZED VIEW public.epitope_11060_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11060_protein_pr;

REFRESH MATERIALIZED VIEW public.epitope_11060_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11060_envelope_pr;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_protein_2k;

REFRESH MATERIALIZED VIEW public.epitope_11060_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11060_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_11060_sfrna1;

REFRESH MATERIALIZED VIEW public.epitope_11060_sfrna2;

REFRESH MATERIALIZED VIEW public.epitope_11060_sfrna3;

REFRESH MATERIALIZED VIEW public.epitope_11060_sfrna4;

-- dengue virus 1
REFRESH MATERIALIZED VIEW public.epitope_11053_polyprotein;

REFRESH MATERIALIZED VIEW public.epitope_11053_anchored_ca;

REFRESH MATERIALIZED VIEW public.epitope_11053_capsid_prot;

REFRESH MATERIALIZED VIEW public.epitope_11053_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11053_protein_pr;

REFRESH MATERIALIZED VIEW public.epitope_11053_membrane_gl;

REFRESH MATERIALIZED VIEW public.epitope_11053_envelope_pr;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_protein_2k;

REFRESH MATERIALIZED VIEW public.epitope_11053_nonstructur;

REFRESH MATERIALIZED VIEW public.epitope_11053_rna_depende;

REFRESH MATERIALIZED VIEW public.epitope_11053_sfrna1;

REFRESH MATERIALIZED VIEW public.epitope_11053_sfrna2;

REFRESH MATERIALIZED VIEW public.epitope_11053_sfrna3;

REFRESH MATERIALIZED VIEW public.epitope_11053_sfrna4;
" | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"


echo "* Update+Refresh of ALL-but-SC2 epitopes completed on database ${database_name} at $(timestamp) ." | tee -a $log_file_path