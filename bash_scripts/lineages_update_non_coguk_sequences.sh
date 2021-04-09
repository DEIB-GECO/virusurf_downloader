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
echo "Check of arguments:
- VIRUSURF DIR: ${virusurf_dir}
- VIRUSURF NORMAL DB: ${database_name};
The program resumes in 10 seconds."
sleep 10
# internal variables
log_file_path=${virusurf_dir}logs/weekly_lineages.log
notifier_script_path=${virusurf_dir}bash_scripts/notify_me.sh


# define utility functions
timestamp() {
  TZ=Europe/Rome date +"%Y_%m_%d__%H_%M_%S" # current time
}
notify_me() { # expects a message as argument
  if test -f "$notifier_script_path"; then
    /bin/bash "$notifier_script_path" "$1"
  else
    echo "virusurf lineage self-update notification not available"  | tee -a $log_file_path
  fi
}
check_exit_code() {    # expects exit status code as argument
  if [ "$1" -ne 0 ]; then
    echo "* Last command terminated abnormally (exit code ${1}). Script interrupted." | tee -a $log_file_path
    notify_me "* $(timestamp) virusurf lineage self-update notification: one of the lineage update steps failed. Check $log_file_path
    for details." | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    exit 1
  fi
}


echo "The output of this script is appended to ${log_file_path}"

# Begin logging of operations
echo "################################          ##################################" | tee -a $log_file_path
echo "####################   Lineage Updater for ViruSurf   @@####################" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic lineage update script started at $(timestamp)" | tee -a $log_file_path



# # Update lineages (of all sequences)
echo "Switching to  conda environment 'vcm'."
source ~/anaconda3_new/etc/profile.d/conda.sh  # finds command conda
check_exit_code "$?"
conda activate vcm
check_exit_code "$?"
cd $virusurf_dir
check_exit_code "$?"
echo "* Begin computing new lineages for ALL SC2 GenBank|NMDC|RefSeq sequences at $(timestamp)" | tee -a $log_file_path
python main.py lineages ${database_name} sars_cov_2 all_except_coguk
check_exit_code "$?"


echo "* Script terminated normally at $(timestamp)" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path

