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
log_file_path=${virusurf_dir}logs/daily_update.log
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
  if [ "$1" -ne 0 ]; then
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

# check if an update was left incomplete or an error occurred in a previous run
if test -f "$ongoing_update_or_error_file_path"; then
    msg="* $(timestamp) virusurf self-update routine: A previous update routine was left incomplete or encountered
    an error. Update interrupted."
    echo "$msg" | tee -a $log_file_path
    notify_me "$msg"
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    echo "" | tee -a $log_file_path
    exit 1
fi

# create update_or_error_flag
touch "$ongoing_update_or_error_file_path"

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



# # Update GISAID
cd $virusurf_dir
echo "* Begin update of GISAID SC2 at $(timestamp)" | tee -a $log_file_path
python main.py import ${database_name_gisaid} gisaid 0 30000
check_exit_code "$?"



# # Update DB with new sequences
cd $virusurf_dir
echo "* Begin update of NCBI SC2 at $(timestamp)" | tee -a $log_file_path
python main.py import ${database_name} sars_cov_2 0 5000
check_exit_code "$?"

# echo "* Begin update of NCBI Bombali at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} bombali_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI MERS at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} mers
# check_exit_code "$?"

# echo "* Begin update of NCBI SC1 at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} sars_cov_1
# check_exit_code "$?"

# echo "* Begin update of NCBI Tai Forest Ebolavirus at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} tai_forest_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI Bundibugyo at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} bundibugyo_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI Reston at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} reston_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI Sudan at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} sudan_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI Zaire at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} zaire_ebolavirus
# check_exit_code "$?"

# echo "* Begin update of NCBI Dengue 4 at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} dengue_4
# check_exit_code "$?"

# echo "* Begin update of NCBI Dengue 3 at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} dengue_3
# check_exit_code "$?"

# echo "* Begin update of NCBI Dengue 2 at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} dengue_2
# check_exit_code "$?"

# echo "* Begin update of NCBI Dengue 1 at $(timestamp)" | tee -a $log_file_path
# python main.py import ${database_name} dengue_1
# check_exit_code "$?"

# # Import COG-UK
echo "* Begin update of COGUK SC2 at $(timestamp)" | tee -a $log_file_path
python main.py import ${database_name} coguk 0 13000
check_exit_code "$?"



# # Refresh materialized views (only 1)
echo "* Refresh of materialized view ${database_name}.nucleotide_variant_annotated at $(timestamp)" | tee -a $log_file_path
psql -U geco -hlocalhost -d ${database_name} -c "REFRESH MATERIALIZED VIEW public.nucleotide_variant_annotated;" | tee -a $log_file_path
check_exit_code "$?"




echo "* Script terminated normally at $(timestamp)" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path
echo "" | tee -a $log_file_path


# remove update_or_error_flag (this happens if the script terminated successfully)
rm "$ongoing_update_or_error_file_path"