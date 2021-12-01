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
log_file_path=${virusurf_dir}logs/daily_update_coguk.log
ongoing_update_or_error_file_path=${virusurf_dir}bash_scripts/ongoing_update_or_error_flag
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
  if [ "$1" -eq 2 ]; then       # code 2 is returned by virusurf_downloader importer when termianted by the user
    echo "* Last python command was terminated by the user (exit code ${1}). Script interrupted." | tee -a $log_file_path
    exit 0
  elif [ "$1" -ne 0 ]; then
    echo "* Last command terminated abnormally (exit code ${1}). Script interrupted." | tee -a $log_file_path
    notify_me "* $(timestamp) virusurf self-update notification: one of the update steps failed. Check $log_file_path
    for details." | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    exit 1
  fi
}


############################################################################

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
echo "######################   Auto Updater for ViruSurf   #######################" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic data update script started at $(timestamp)" | tee -a $log_file_path



# # Import COG-UK
cd $virusurf_dir
echo "* Begin update of COGUK SC2 at $(timestamp)" | tee -a $log_file_path
python main.py import ${database_name} coguk 0 60000
check_exit_code "$?"


# # Refresh materialized views (only 1)
echo "* Refresh of materialized view ${database_name}.nucleotide_variant_annotated at $(timestamp)" | tee -a $log_file_path
psql -U geco -hlocalhost -d ${database_name} -c "REFRESH MATERIALIZED VIEW public.nucleotide_variant_annotated;" | tee -a $log_file_path
check_exit_code "$?"



echo "* Script terminated normally at $(timestamp)" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path
