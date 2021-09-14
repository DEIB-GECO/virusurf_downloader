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

# define utility functions
timestamp() {
  TZ=Europe/Rome date +"%Y_%m_%d__%H_%M_%S" # current time
}
check_exit_code() {    # expects exit status code as argument
  if [ "$1" -ne 0 ]; then
    echo "* Last command terminated abnormally (exit code ${1}). Script interrupted."
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
cd $virusurf_dir
echo "Stopping import of genbank sars cov 2 on ${database_name}"
python main.py stop import ${database_name} sars_cov_2
echo "Stopping import of coguk sars cov 2 on ${database_name}"
python main.py stop import ${database_name} coguk
echo "Stopping import of gisaid on ${database_name_gisaid}"
python main.py stop import ${database_name_gisaid} gisaid

