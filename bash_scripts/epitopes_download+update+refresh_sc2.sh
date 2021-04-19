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

log_file_path=${virusurf_dir}logs/epitope_sc2_updates.log
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
echo "##################   Epitope SC2 Updater for ViruSurf   ####################" | tee -a $log_file_path
echo "#######                - downloads updates from IEDB                 #######" | tee -a $log_file_path
echo "#######                - deletes outdated SC2 epitopes               #######" | tee -a $log_file_path
echo "#######                - import SC2-related epitopes                 #######" | tee -a $log_file_path
echo "#######                - refresh SC2-related views for EpiSurf       #######" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic epitopes update script started at $(timestamp)" | tee -a $log_file_path


# Update epitope files from IEDB
if [ "$enable_download" == "yes" ] || [ "$enable_download" == "true" ]; then
  echo "* Download of updated source files from IEDB at $(timestamp) ... " | tee -a $log_file_path
  cd $virusurf_dir
  python main.py download epitopes
  check_exit_code "$?"
fi


# Delete old epitopes
echo "* Deletion of previous SC2 epitopes at $(timestamp) ..." | tee -a $log_file_path
cd $virusurf_dir
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -c "
-- delete epitope fragments of this virus
WITH which_virus AS (
    SELECT virus_id
    FROM virus
    WHERE taxon_id = 2697049 -- SC2 taxon_id
)
DELETE FROM epitope_fragment AS ef
USING epitope AS e
WHERE ef.epitope_id = e.epitope_id AND e.virus_id = (SELECT which_virus.virus_id FROM which_virus);

-- delete epitopes of this virus
WITH which_virus AS (
    SELECT virus_id
    FROM virus
    WHERE taxon_id = 2697049 -- SC2 taxon_id
)
DELETE FROM epitope AS e
WHERE e.virus_id = (SELECT which_virus.virus_id FROM which_virus);
" | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"

# Insert SC2 epitopes
echo "* Begin import of epitopes for virus sars_cov_2 at $(timestamp)" | tee -a $log_file_path
cd $virusurf_dir
python main.py epitopes ${database_name} sars_cov_2
check_exit_code "$?"


# Refresh materialized views for episurf
echo "* Refresh SC2 eptiope variants and info materialized view at $(timestamp)"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -c "
REFRESH MATERIALIZED VIEW public.epitope_2697049_orf1ab_poly;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp12_rna_d;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp13_helic;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp14_3_to_;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp15_endor;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp16_2_o_r;

REFRESH MATERIALIZED VIEW public.epitope_2697049_orf1a_polyp;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp1_leader;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp2;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp3;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp4;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp5_3c_lik;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp6;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp7;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp8;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp9;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp10;

REFRESH MATERIALIZED VIEW public.epitope_2697049_nsp11;

REFRESH MATERIALIZED VIEW public.epitope_2697049_spike_surfa;

REFRESH MATERIALIZED VIEW public.epitope_2697049_ns3_orf3a_p;

REFRESH MATERIALIZED VIEW public.epitope_2697049_e_envelope_;

REFRESH MATERIALIZED VIEW public.epitope_2697049_m_membrane_;

REFRESH MATERIALIZED VIEW public.epitope_2697049_ns6_orf6_pr;

REFRESH MATERIALIZED VIEW public.epitope_2697049_ns7a_orf7a_;

REFRESH MATERIALIZED VIEW public.epitope_2697049_ns7b_orf7b;

REFRESH MATERIALIZED VIEW public.epitope_2697049_ns8_orf8_pr;

REFRESH MATERIALIZED VIEW public.epitope_2697049_n_nucleocap;

REFRESH MATERIALIZED VIEW public.epitope_2697049_orf10_prote;
" | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"


echo "* Update of SC2 epitopes completed on database ${database_name} at $(timestamp)." | tee -a $log_file_path