#!/bin/bash

# make sure this script is executable:  chmod 755 myscript.sh
# set your user as owner of the script: chown myusername: myscript.sh

echo "Assumptions for the proper execution of this script:
- command line arguments: <virusurf_directory_path> <virusurf_updating_databases_path>
- the postgres password file contains one entry for each of the given databases formatted as host:port:dat_name:user:psw;
- a conda environment named 'vcm' exists and it has the packages necessary for virusurf_downlaoder already installed;"
sleep 10


# input parameters
virusurf_dir=${1}/
virusurf_updating_databases_path=$2
database_name=$(head -1 "$virusurf_updating_databases_path")
database_name_gisaid=$(head -2 "$virusurf_updating_databases_path" | tail -1)
echo "Check of arguments:
- VIRUSURF DIR: ${virusurf_dir}
- VIRUSURF NORMAL DB: ${database_name};
- VIRUSURF GISAID DB ${database_name_gisaid};
The program resumes in 10 seconds."
sleep 10
# internal variables
log_file_path=${virusurf_dir}logs/weekly_overlaps.log
ongoing_update_or_error_file_path=${virusurf_dir}bash_scripts/ongoing_update_or_error_flag
ongoing_overlap_or_error_file_path=${virusurf_dir}bash_scripts/ongoing_overlap_or_error_flag
notifier_script_path=${virusurf_dir}bash_scripts/notify_me.sh
db_overlaps_configuration_path="${virusurf_dir}db_config/db_overlaps_configuration.json"
db_overlaps_configuration_temp_path="${virusurf_dir}db_config/db_overlaps_configuration_temp.json"


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
    notify_me "* $(timestamp) virusurf self-update notification: one of the overlps steps failed. Check $log_file_path
    for details. RESET OF db_overlaps_configuration.json REQUIRED (moved to db_overlaps_configuration_temp.json." | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    mv "$db_overlaps_configuration_path" "$db_overlaps_configuration_temp_path"
    exit
  fi
}


############################################################################

# check if an update was left incomplete or an error occurred in a previous run
if test -f "$ongoing_update_or_error_file_path" || test -f "$ongoing_overlap_or_error_file_path"; then
    msg="* $(timestamp) virusurf self-update routine: A previous update/overlap routine was left incomplete or encountered
    an error. Update interrupted."
    echo "$msg" | tee -a $log_file_path
    notify_me "$msg"
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    exit
fi

# create update_or_error_flag
touch "$ongoing_overlap_or_error_file_path"

# switch to conda environment vcm
echo "Switching to  conda environment 'vcm'."
source ~/miniconda3/etc/profile.d/conda.sh  # finds command conda
check_exit_code "$?"
conda activate vcm
check_exit_code "$?"

##############################################################################


echo "The output of this script is appended to ${log_file_path}"

# Begin logging of operations
echo "################################          ##################################" | tee -a $log_file_path
echo "######################   Auto Updater for ViruSurf   #######################" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic overlaps-checking script started at $(timestamp)" | tee -a $log_file_path


# update db_overlaps_config JSON
cd $virusurf_dir
echo "* Setup of db_overlaps_config JSON at $(timestamp)" | tee -a $log_file_path
python find_n_replace_in_file.py "$db_overlaps_configuration_path" "vcm_gisaid_du" "${database_name_gisaid}"
check_exit_code "$?"
python find_n_replace_in_file.py "$db_overlaps_configuration_path" "vcm_du" "${database_name}"
check_exit_code "$?"


# # Update DB with new overlaps
cd $virusurf_dir
echo "* Begin overlaps GenBank-GISAID at $(timestamp)" | tee -a $log_file_path
python main.py overlaps genbank_gisaid True
check_exit_code "$?"


echo "* Begin overlaps COGUK-GISAID at $(timestamp)" | tee -a $log_file_path
python main.py overlaps coguk_gisaid True
check_exit_code "$?"

echo "* Begin overlaps NMDC-GISAID at $(timestamp)" | tee -a $log_file_path
python main.py overlaps gisaid_nmdc True
check_exit_code "$?"



echo "* Script terminated normally at $(timestamp)" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path



# restore db_overlaps_config JSON
cd $virusurf_dir
echo "* Restore of db_overlaps_config JSON at $(timestamp)" | tee -a $log_file_path
python find_n_replace_in_file.py "$db_overlaps_configuration_path" "${database_name_gisaid}" "vcm_gisaid_du"
check_exit_code "$?"
python find_n_replace_in_file.py "$db_overlaps_configuration_path" "${database_name}" "vcm_du"
check_exit_code "$?"


# remove update_or_error_flag (this happens if the script terminated successfully)
rm "$ongoing_overlap_or_error_file_path"