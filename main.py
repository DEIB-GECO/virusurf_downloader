import os
import signal
import sys
from loguru import logger

import logger_settings
from data_sources.ncbi_any_virus.ncbi_importer import import_samples_into_vcm
from data_sources.ncbi_any_virus.settings import known_settings
import data_sources.nmdc.procedure as nmdc
import stats_module
from logger_settings import setup_logger
from db_config import read_db_import_configuration as import_config
from os.path import abspath
from generate_fasta import generate_fasta
from datetime import date
from locations import get_local_folder_for, FileType
import psutil
from datetime import datetime, timedelta

log_file_keyword = ""

#   #################################       PROGRAM ARGUMENTS   ##########################
ncbi_virus_names = sorted(known_settings.keys())
wrong_arguments_message = "Accepted parameters:\n" \
                          "\timport <db_name> coguk|gisaid|" + "|".join(ncbi_virus_names) + "\n" \
                          "\tdownload epitopes" + "\n" \
                          "\tepitopes <db_name> " + "|".join(ncbi_virus_names) + "\n" \
                          "\toverlaps <source1>_<source2> <commit to the DB ?>\n" \
                          "\tlineages <db_name> " + "|".join(ncbi_virus_names) + "all|only_null_lineages|all_except_coguk" + "\n" \
                          "\tstop <import|download|epitopes|overlaps|lineages> <optional_arg> <optional_arg>"
wrong_arguments_message += 'When action is "import", the source name can be optionally followed by a range of samples to import as <min> (included) <max> (excluded).\n'


def cancel_an_old_run_with_parameters(args: list):
    logger.info(f"Searching an older instance running with keywords: python + {' '.join(args)}")
    process_command_words = set(args)
    result = []
    for process in psutil.process_iter():
        try:
            if 'python' in process.name() and process_command_words.issubset(set(process.cmdline())):
                result.append((process.pid, process.cmdline(), process.ppid()))
        except psutil.AccessDenied:
            continue

    ## terminate old run
    previous_instances = [x for x in result if x[0] != os.getpid()]
    if len(previous_instances) > 1:
        # the main process can spawn many worker processes... -> remove children processes
        parent_pids = {y[2] for y in previous_instances}
        previous_instances = [x for x in previous_instances if x[0] in parent_pids]
    if previous_instances:
        if len(previous_instances) > 1:
            raise SystemError(f"More than one old instance is still running. Maybe there is a bug. "
                              f"Current time is {datetime.now()} UTC.")
        else:
            process_to_terminate = psutil.Process(previous_instances[0][0])
            logger.info(f"Terminating old instance: {' '.join(process_to_terminate.cmdline())} with PID "
                        f"{process_to_terminate.pid}")
            process_to_terminate.send_signal(signal.SIGINT)
            SECONDS_TO_WAIT_FOR_TERMINATION = 60 * 60 * 4
            time_to_wait_timedelta_string = timedelta(seconds=SECONDS_TO_WAIT_FOR_TERMINATION)
            try:
                logger.info(f"Waiting at most {time_to_wait_timedelta_string} (h:m:s) for process PID "
                            f"{process_to_terminate.pid} to terminate itself.")
                process_to_terminate.wait(SECONDS_TO_WAIT_FOR_TERMINATION)     # Keyboard Interrupt handled outside
            except psutil.TimeoutExpired:
                raise TimeoutError(f"Process {' '.join(process_to_terminate.cmdline())} with PID "
                                   f"{process_to_terminate.pid} is still alive after "
                                   f"{time_to_wait_timedelta_string} (h:m:s). Maybe there is a bug. "
                                   f"Current time is {datetime.now()} LOCAL TIME.")
    else:
        logger.info(f"No old instance found.")


def start():
    try:
        action = sys.argv[1].lower()
        if action in ['import', 'epitopes', 'lineages']:  # get database name if action is not 'overlaps' or 'download'
            db_name = sys.argv[2]
            import_config.set_db_name(db_name)
        if action == 'import':
            source = sys.argv[3].lower()
            log_file_keyword = source
            try:
                _from = int(sys.argv[4])
                to = int(sys.argv[5])
                if not _from < to:
                    logger.error('Optional parameter <min> must be less than <max>')
                    sys.exit(1)
            except IndexError:
                _from = None
                to = None
        elif 'epitopes' in action:
            _epitope_target = sys.argv[3]
            known_settings[_epitope_target]  # raises KeyError if the source is not recognized
            log_file_keyword = f"epi_{_epitope_target}"
        elif 'lineages' in action:
            _fasta_target = sys.argv[3]
            known_settings[_fasta_target]  # raises KeyError if the source is not recognized
            _method = sys.argv[4]
            if _method == 'all':
                _only_null_lineages = False
                _data_source = None
            elif _method == 'only_null_lineages':
                _only_null_lineages = True
                _data_source = None
            elif _method == 'all_except_coguk':
                _only_null_lineages = False
                _data_source = ('GenBank', 'NMDC', 'RefSeq')
            else:
                raise ValueError('last argument must be one of all|only_null_lineages|all_except_coguk')
            log_file_keyword = f"lineages_{_fasta_target}"
        elif 'overlaps' in action:
            _overlap_target = sys.argv[2]
            log_file_keyword = f"overlaps_{_overlap_target}"
        elif 'download' in action:
            _what_to_download = sys.argv[2]
            log_file_keyword = f"download_{_what_to_download}"
        elif 'stop' in action:
            _commands_of_run_to_be_stopped = sys.argv[2:5]
            log_file_keyword = f"stop_{'_'.join(_commands_of_run_to_be_stopped)}"
        else:
            log_file_keyword = action
    except (IndexError, KeyError):
        logger.error(wrong_arguments_message)
        sys.exit(1)

    #   ###################################      SETUP LOGGER    ##############################
    setup_logger(log_file_keyword)

    #   #############################     CANCEL PREVIOUS RUN (IF ANY)      ###################
    # find old runs by looking to processes executing a command like this one, up to the 4th argument
    # for example: main.py import db_name sars_cov_2
    try:
        cancel_an_old_run_with_parameters(sys.argv[:4])
    except TimeoutError as e:
        error_msg = f"FATAL ERROR.\nTimeoutError {e}"
        logger.error(error_msg)
        logger_settings.send_message(error_msg, block=True)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Operation aborted. Launch of a new instance canceled.")
        return

    #   ###################################     PERFORM <action>       ########################
    try:
        if 'download' in action:
            _what_to_download = sys.argv[2]
            if 'epitopes' in _what_to_download:
                from epitopes import download_epitope_data
                download_epitope_data()
        elif 'epitopes' in action:
            from epitopes import import_epitopes
            # noinspection PyUnboundLocalVariable
            virus_import_parameters = known_settings[_epitope_target]
            virus_txid = virus_import_parameters["virus_taxon_id"]
            import_epitopes(virus_txid)
        elif 'import' in action:
            # noinspection PyUnboundLocalVariable
            if source in ['coguk', 'cog-uk']:
                from data_sources.coguk_sars_cov_2.procedure import run as run_coguk
                run_coguk(from_sample=_from, to_sample=to)
            elif source == 'gisaid':
                from data_sources.gisaid_sars_cov_2.procedure import run as run_gisaid
                run_gisaid(from_sample=_from, to_sample=to)
            elif source in known_settings.keys():
                import_samples_into_vcm(source, from_sample=_from, to_sample=to)
            elif 'nmdc' in source:
                nmdc.import_samples_into_vcm()
            else:
                logger.error(f'the argument {source} is not recognised.\n' + wrong_arguments_message)
        elif 'lineages' in action:
            from lineages import update_pangolin, call_pangolin, update_db_with_computed_lineages
            # create fasta with the nucleotide sequences to analyze with pangolin
            virus_import_parameters = known_settings[_fasta_target]
            virus_txid = virus_import_parameters["virus_taxon_id"]
            virus_folder = virus_import_parameters["generated_dir_name"]
            fasta_name = f'{_fasta_target}_{_method}_{date.today().strftime("%Y-%b-%d")}.fasta'
            fasta_path = generate_fasta(virus_txid, virus_folder, fasta_name, _only_null_lineages, _data_source)
            # use pangoling to compute lineages
            update_pangolin()
            pangolin_output_path = abspath(get_local_folder_for(virus_folder, FileType.Fasta)
                                           + 'lineages_pangolin_output.csv')
            call_pangolin(fasta_path, pangolin_output_path)
            # load lineages into db
            update_db_with_computed_lineages(pangolin_output_path)
        elif 'overlaps' in action:
            from overlaps import overlaps_controller
            overlaps_controller.run()
        elif 'stop' in action:
            _commands_of_run_to_be_stopped.insert(0, sys.argv[0])  # add this module's name
            cancel_an_old_run_with_parameters(_commands_of_run_to_be_stopped)
        else:
            logger.error(f'the argument {action} is not recognised.\n' + wrong_arguments_message)
    except:
        logger.exception(
            'FATAL ERROR')  # this is just to make sure the exception is written to the log file before crashing
        sys.exit(1)
    finally:
        if 'import' in action:
            stats_module.check_samples_imported()
        logger.complete()


if __name__ == '__main__':
    start()
