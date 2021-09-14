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
echo "filling epitope tables for virus dengue_1"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_dengue_1.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus dengue_2"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_dengue_2.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus dengue_3"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_dengue_3.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus dengue_4"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_dengue_4.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus mers"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_mers.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus reston_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_reston_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus sars_cov_1"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_sars_cov_1.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus sudan_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_sudan_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus tai_forest_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_tai_forest_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus zaire_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_zaire_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus bombali_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_bombali_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"
echo "filling epitope tables for virus bundibugyo_ebolavirus"
psql -v ON_ERROR_STOP=1 -h localhost -U geco -d ${database_name} -f ${virusurf_dir}sql_scripts/epitope_tables/refresh_table/refresh_bundibugyo_ebolavirus.sql | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"


echo "* Update+Refresh of ALL-but-SC2 epitopes completed on database ${database_name} at $(timestamp) ." | tee -a $log_file_path