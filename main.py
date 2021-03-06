import sys
from loguru import logger
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

log_file_keyword = ""

#   #################################       PROGRAM ARGUMENTS   ##########################
ncbi_virus_names = sorted(known_settings.keys())
wrong_arguments_message = "Accepted parameters:\n" \
                          "\timport <db_name> coguk|gisaid|" + "|".join(ncbi_virus_names) + "\n" \
                          "\tdownload epitopes" + "\n" \
                          "\tepitopes <db_name> " + "|".join(ncbi_virus_names) + "\n" \
                          "\toverlaps <source1>_<source2> <commit to the DB ?>\n" \
                          "\tlineages <db_name> " + "|".join(ncbi_virus_names) + "all|only_null_lineages|all_except_coguk" + "\n"
wrong_arguments_message += 'When action is "import", the source name can be optionally followed by a range of samples to import as <min> (included) <max> (excluded).\n'

try:
    action = sys.argv[1].lower()
    if action in ['import', 'epitopes', 'lineages']:     # get database name if action is not 'overlaps' or 'download'
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
        known_settings[_epitope_target]     # raises KeyError if the source is not recognized
        log_file_keyword = f"epi_{_epitope_target}"
    elif 'lineages' in action:
        _fasta_target = sys.argv[3]
        known_settings[_fasta_target]       # raises KeyError if the source is not recognized
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
    else:
        log_file_keyword = action
except (IndexError, KeyError):
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
setup_logger(log_file_keyword)


#   ###################################     PERFORM <action>       ###############
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
            logger.error(f'the argument {source} is not recognised.\n'+wrong_arguments_message)
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
    else:
        logger.error(f'the argument {action} is not recognised.\n' + wrong_arguments_message)
except:
    logger.exception('FATAL ERROR')  # this is just to make sure the exception is written to the log file before crashing
    sys.exit(1)
finally:
    if 'import' in action:
        stats_module.check_samples_imported()
    logger.complete()
