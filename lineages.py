from datetime import date
from subprocess import call, check_output, CalledProcessError
import subprocess
from loguru import logger
from locations import remove_file
from os.path import sep, abspath
from db_config.read_db_import_configuration import get_database_config_params


conda_install_path = '/home/metadata/anaconda3_new/etc/profile.d/conda.sh'
logger.warning(f'Assume to find conda installed at path: {conda_install_path}. If not, change the address in module '
               f'lineages.py')
logger.warning('This action requires pangolin to be installed in a pangolin conda environment.')


def update_pangolin():
    shell_cmds = [
        f'source {conda_install_path}',
        'conda activate pangolin',
        'pangolin --update'
    ]

    shell_cmds = ' && '.join(shell_cmds)
    logger.trace(f"executing bash commands {shell_cmds}")
    try:
        output = check_output(shell_cmds, shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)  # exceptions are caught externally
        logger.info(f"Output of pangolin --update:\n{output.decode('utf-8')}")
    except OSError as e:
        logger.error(f"the process updating pangolin raised an exception")
        raise e
    except CalledProcessError as e:
        ret_code = e.returncode
        output = e.output
        logger.error(f'the process updating pangolin returned with non-zero exit code {ret_code} and output:\n {output}')
        raise e


def call_pangolin(input_fasta_path: str, pangolin_output_path: str):
    pangolin_output_dir = pangolin_output_path[:pangolin_output_path.rindex(f'{sep}')]
    pangolin_output_file_name = pangolin_output_path[pangolin_output_path.rindex(f'{sep}')+1:]
    shell_cmds = [
        f'source {conda_install_path}',
        'conda activate pangolin',
        f'pangolin {input_fasta_path} -o {pangolin_output_dir} --outfile {pangolin_output_file_name}'
    ]

    shell_cmds = ' && '.join(shell_cmds)
    logger.trace(f"executing bash command {shell_cmds} ")
    try:
        ret_code = call(shell_cmds, shell=True, executable='/bin/bash')  # exceptions are caught externally
    except OSError as e:
        logger.error(f"the process running pangolin raised an exception")
        raise e
    if ret_code < 0:
        raise ChildProcessError(f"the process running pangolin was terminated by signal {-ret_code}")
    elif ret_code != 0:
        raise ChildProcessError(f"the process running pangolin returned with non-zero exit code ({ret_code})")


def update_db_with_computed_lineages(path_to_pangolin_output: str):
    sql_loader_script = abspath(f'.{sep}sql_scripts{sep}load_lineages.sql')

    # create a temporary working sql script
    tmp_script_path = abspath(f'.{sep}sql_scripts{sep}load_lineages_{date.today().strftime("%Y-%b-%d")}.sql')
    with open(tmp_script_path, mode='w') as temp_script:
        with open(sql_loader_script, mode='r') as template_script:
            for line in template_script.readlines():
                new_line = line.replace('pathtocsvfile', path_to_pangolin_output)
                temp_script.write(new_line)

    # execute the sql script with psql
    db_params = get_database_config_params()
    db_user = db_params["db_user"]
    db_name = db_params["db_name"]
    psql_cmd = f'psql -v ON_ERROR_STOP=1 -h localhost -U {db_user} -d {db_name} -f {tmp_script_path}'
    logger.trace(f"executing bash command {psql_cmd}")
    try:
        output = check_output(psql_cmd, shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
        logger.info(f"Output of {psql_cmd}:\n{output.decode('utf-8')}")
        # remove the temporary working copy
        remove_file(tmp_script_path)
    except OSError as e:
        logger.error(f"the psql process updating lineages raised an exception while executing {tmp_script_path}")
        raise e
    except CalledProcessError as e:
        ret_code = e.returncode
        output = e.output
        logger.error(f'the psql process updating lineages returned with non-zero exit code {ret_code}  while '
                     f'executing {tmp_script_path}. Output is:\n {output}')
        raise e
