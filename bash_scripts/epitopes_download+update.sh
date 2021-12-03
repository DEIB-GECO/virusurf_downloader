#!/bin/bash

# make sure this script is executable:  chmod 755 myscript.sh
# set your user as owner of the script: chown myusername: myscript.sh

echo "Assumptions for the proper execution of this script:
- command line arguments: <virusurf_directory_path> <target_database_name>;
- the postgres password file contains one entry for the given database formatted as host:port:dat_name:user:psw;
- a conda environment named 'vcm' exists and it has the packages necessary for virusurf_downlaoder already installed;"
echo "THIS SCRIPT DOES NOT UPDATE EPITOPE ASSOCIATED MATERIALIZED VIEWS"

# input parameters
virusurf_dir=${1}/
database_name=$2
echo "Check of arguments:
- VIRUSURF DIR: ${virusurf_dir};
- VIRUSURF DB: ${database_name};
The program resumes in 10 seconds."
sleep 10

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
    exit 1
  fi
}

echo "The output of this script is appended to ${log_file_path}"


# switch to conda environment vcm
echo "Switching to  conda environment 'vcm'."
source ~/anaconda3_new/etc/profile.d/conda.sh  # finds command conda
check_exit_code "$?"
conda activate vcm
check_exit_code "$?"

# Begin logging of operations
echo "################################          ##################################" | tee -a $log_file_path
echo "####################    Epitope Updater for ViruSurf   #####################" | tee -a $log_file_path
echo "################   (does NOT affect views for EpiSurf)   ###################" | tee -a $log_file_path
echo "################################          ##################################" | tee -a $log_file_path
echo "* ViruSurf automatic epitopes update script started at $(timestamp)" | tee -a $log_file_path



# Update epitope files from IEDB
echo "* Download of updated source files from IEDB..." | tee -a $log_file_path
cd $virusurf_dir
python main.py download epitopes
check_exit_code "$?"

# Delete old epitopes
echo "* Deletion of previous epitopes..." | tee -a $log_file_path
cd $virusurf_dir
psql -h localhost -U geco -d ${database_name} -c "TRUNCATE TABLE epitope_fragment, epitope;" | tee -a $log_file_path
check_exit_code "${PIPESTATUS[0]}"

# Insert new epitopes
cd $virusurf_dir
echo "* Begin import of epitopes for virus bundibugyo_ebolavirus" | tee -a $log_file_path
python main.py epitopes ${database_name} bundibugyo_ebolavirus
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_1" | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_1
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_2" | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_2
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_3" | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_3
check_exit_code "$?"

echo "* Begin import of epitopes for virus dengue_4" | tee -a $log_file_path
python main.py epitopes ${database_name} dengue_4
check_exit_code "$?"

echo "* Begin import of epitopes for virus mers" | tee -a $log_file_path
python main.py epitopes ${database_name} mers
check_exit_code "$?"

echo "* Begin import of epitopes for virus sars_cov_2" | tee -a $log_file_path
python main.py epitopes ${database_name} sars_cov_2
check_exit_code "$?"

echo "* Begin import of epitopes for virus sars_cov_1" | tee -a $log_file_path
python main.py epitopes ${database_name} sars_cov_1
check_exit_code "$?"

echo "* Begin import of epitopes for virus sudan_ebolavirus" | tee -a $log_file_path
python main.py epitopes ${database_name} sudan_ebolavirus
check_exit_code "$?"

echo "* Begin import of epitopes for virus zaire_ebolavirus" | tee -a $log_file_path
python main.py epitopes ${database_name} zaire_ebolavirus
check_exit_code "$?"

echo "* Import of epitopes completed on database ${database_name}." | tee -a $log_file_path
echo "THIS SCRIPT DID NOT UPDATE EPITOPE ASSOCIATED MATERIALIZED VIEWS" | tee -a $log_file_path


