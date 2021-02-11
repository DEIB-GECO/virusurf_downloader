#!/bin/bash

# make sure this script is executable:  chmod 755 myscript.sh
# set your user as owner of the script: chown myusername: myscript.sh

echo "This script expects the database_name as parameter."
database_name=$1
echo "Epitopes are being updated on database ${database_name}.
Make sure that the password for database ${database_name} is present in file .pgpass as host:port:dat_name:user:psw"
sleep 10


user_home_path=/home/alfonsi/
virusurf_dir=${user_home_path}virusurf_downloader/
vcm_db_backups_dir=${user_home_path}ViruSurf_backups/
log_file_path=${virusurf_dir}logs/epitope_updates.log

# define utility functions
timestamp() {
  TZ=Europe/Rome date +"%Y_%m_%d__%H_%M_%S" # current time
}
check_exit_code() {    # expects exit status code as argument
  if [ "$1" -ne 0 ]; then
    echo "* Last command terminated abnormally (exit code ${1}). Script interrupted." | tee -a $log_file_path
    echo ""
    echo ""
    echo ""
    exit
  fi
}
backup_db() { # expect the full_path of the backup to be created.
  echo "* Backup of DB overwriting file $1. PSQL OUTPUT BELOW" | tee -a $log_file_path
  pg_dump -U geco -hlocalhost ${database_name} > $1
}

echo "The output of this script is appended to ${log_file_path}"

# Begin logging of operations
#echo "################################          ##################################" | tee -a $log_file_path
#echo "######################   Auto Updater for ViruSurf   #######################" | tee -a $log_file_path
#echo "################################          ##################################" | tee -a $log_file_path
#echo "* ViruSurf automatic data update script started at $(timestamp)" | tee -a $log_file_path


# Do a backup first
#echo "* Backup of database ${database_name} in progress..." | tee -a $log_file_path
#cd $vcm_db_backups_dir
#backup_db "${vcm_db_backups_dir}${database_name}_@_$(timestamp).dump"
#check_exit_code "$?"
#
## Delete old epitopes
#echo "* Deletion of previous epitopes..." | tee -a $log_file_path
#psql -h localhost -U geco -d ${database_name} -c "TRUNCATE TABLE epitope, epitope_fragment;"
#check_exit_code "$?"
#
## Update epitope files from IEDB
#echo "* Download of updated source files from IEDB..." | tee -a $log_file_path
#cd $virusurf_dir
#python main.py download epitopes
#check_exit_code "$?"

# Insert new epitopes
cd $virusurf_dir
#echo "* Begin import of epitopes for virus bundibugyo_ebolavirus" | tee -a $log_file_path
#python main.py epitopes ${database_name} bundibugyo_ebolavirus
#check_exit_code "$?"
#
#echo "* Begin import of epitopes for virus dengue_1" | tee -a $log_file_path
#python main.py epitopes ${database_name} dengue_1
#check_exit_code "$?"
#
#echo "* Begin import of epitopes for virus dengue_2" | tee -a $log_file_path
#python main.py epitopes ${database_name} dengue_2
#check_exit_code "$?"
#
#echo "* Begin import of epitopes for virus dengue_3" | tee -a $log_file_path
#python main.py epitopes ${database_name} dengue_3
#check_exit_code "$?"
#
#echo "* Begin import of epitopes for virus dengue_4" | tee -a $log_file_path
#python main.py epitopes ${database_name} dengue_4
#check_exit_code "$?"

echo "* Begin import of epitopes for virus mers" | tee -a $log_file_path
python main.py epitopes ${database_name} mers
check_exit_code "$?"

echo "* Begin import of epitopes for virus sars_cov_2" | tee -a $log_file_path
python main.py epitopes ${database_name} sars_cov_2
check_exit_code "$?"

echo "* Begin import of epitopes for virus sars_cov_1" | tee -a $log_file_path
python main.py epitopes ${database_name} sars_cov_2
check_exit_code "$?"

echo "* Begin import of epitopes for virus sudan_ebolavirus" | tee -a $log_file_path
python main.py epitopes ${database_name} sudan_ebolavirus
check_exit_code "$?"

echo "* Begin import of epitopes for virus zaire_ebolavirus" | tee -a $log_file_path
python main.py epitopes ${database_name} zaire_ebolavirus
check_exit_code "$?"

echo "* Update of epitopes completed on database ${database_name}." | tee -a $log_file_path

