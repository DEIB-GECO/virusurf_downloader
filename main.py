import sys
from loguru import logger
import database_tom
from data_sources.ncbi_any_virus.ncbi_importer import prepared_parameters, import_samples_into_vcm
import data_sources.nmdc.procedure as nmdc
from Bio import Entrez
from tqdm import tqdm
import stats_module
import warnings
import notifications_module

Entrez.email = "Your.Name.Here@example.org"
log_file_keyword = ""


#   #################################       PROGRAM ARGUMENTS   ##########################
wrong_arguments_message = 'The module main.py always expects at least these arguments: ' \
                          'db_name, recreate_db?, db_user, db_password, db_port, <action>\n' \
                          'Acceptable values for action are\n:' \
                          '\timport <name of the source> <optional range of sequences (e.g. 0 100) for ncbi and coguk>\n' \
                          '\tepitopes <name of the source>\n' \
                          '\toverlaps <source1>_<source2> <commit to the DB ?>\n' \
                          '\tlineages <lineages TSV file path>' \
                          '\tindexes\n' \
                          '\tviews\n' \
                          '\tchimera_sequences\n' \
                          '\tAcceptable source names are:\n' \
                          '\tcoguk\n' \
                          '\tgisaid\n'
for p in prepared_parameters.keys():
    wrong_arguments_message += f'\t{p}\n'
wrong_arguments_message += 'When action is "import", the source name can be optionally followed by a range of samples to import as <min> (included) <max> (excluded).\n'
# noinspection PyBroadException
try:
    db_name = sys.argv[1]
    _db_recreate = sys.argv[2].lower()
    db_recreate = True if _db_recreate == 'yes' or _db_recreate == 'true' else False
    db_user = sys.argv[3]
    db_password = sys.argv[4]
    db_port = sys.argv[5]
    action = sys.argv[6].lower()
    if action == 'import':
        source = sys.argv[7].lower()
        log_file_keyword = source
        try:
            _from = int(sys.argv[8])
            to = int(sys.argv[9])
            if not _from < to:
                logger.error('Optional parameter <min> must be less than <max>')
                sys.exit(1)
        except Exception:
            _from = None
            to = None
    elif 'epitope' in action:
        _epitope_target = sys.argv[7]
        log_file_keyword = f"epi_{_epitope_target}"
    elif 'lineage' in action:
        _lineages_tsv_file_path = sys.argv[7]
        log_file_keyword = f"lineages"
    elif 'fasta' in action:
        _fasta_target = sys.argv[7]
        log_file_keyword = f"fasta_{_fasta_target}"
    elif 'overlap' in action:
        _overlap_target = sys.argv[7]
        log_file_keyword = f"overlaps_{_overlap_target}"
    else:
        log_file_keyword = action
except Exception:
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
logger.remove()  # removes default logger to stderr with level DEBUG
# on console print from level INFO on
logger.add(sink=lambda msg: tqdm.write(msg, end=''),
           level='TRACE',
           format="<green>{time:YYYY-MM-DD HH:mm:ss Z}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True,
           enqueue=True)
# log to file any message of any security level
logger.add("./logs/log_"+log_file_keyword+"_{time}.log",
           level='TRACE',
           rotation='100 MB',
           format="<green>{time:YYYY-MM-DD HH:mm:ss Z}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=False,
           backtrace=True,
           diagnose=True,
           enqueue=True)
notifications_module.setup_any_additional_error_notifiers()
# redirect warnings
def customwarn(message, category, filename, lineno, file=None, line=None):
    logger.warning(warnings.formatwarning(message, category, filename, lineno))
warnings.showwarning = customwarn

logger.info(f"main.py {' '.join(sys.argv[1:])}")

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
# init database
if 'overlaps' not in action:
    database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=db_recreate)

#   ###################################     MAIN OPERATION       ###############
try:
    if 'index' in action:
        database_tom.create_indexes()
    elif 'view' in action:
        database_tom.create_views()
    elif 'chimera_sequence' in action:
        database_tom.disambiguate_chimera_sequences()
    elif 'epitope' in action:
        from epitopes import import_epitopes
        if 'bombali' in _epitope_target or 'reston' in _epitope_target or 'forest' in _epitope_target:
            raise ValueError(f'No epitopes available for virus {_epitope_target}.')
        virus_import_parameters = prepared_parameters.get(_epitope_target)
        if not virus_import_parameters:
            raise ValueError(f'{_epitope_target} is not recognised as an importable virus')
        virus_txid = virus_import_parameters[1]
        import_epitopes(virus_txid)
    elif 'lineage' in action:
        from load_lineages import load
        load(_lineages_tsv_file_path)
    elif 'import' in action:
        if source in ['coguk', 'cog-uk']:
            from data_sources.coguk_sars_cov_2.procedure import run as run_coguk
            run_coguk(from_sample=_from, to_sample=to)
        elif source == 'gisaid':
            from data_sources.gisaid_sars_cov_2.procedure import run as run_gisaid
            run_gisaid(from_sample=_from, to_sample=to)
        elif source in prepared_parameters.keys():
            import_samples_into_vcm(*prepared_parameters[source], from_sample=_from, to_sample=to)
        elif 'nmdc' in source:
            nmdc.import_samples_into_vcm()
        else:
            logger.error(f'the argument {source} is not recognised.\n'+wrong_arguments_message)
    elif 'fasta' in action:
        from generate_fasta import generate_fasta
        virus_import_parameters = prepared_parameters.get(_fasta_target)
        if not virus_import_parameters:
            raise ValueError(f'{_fasta_target} is not recognised as an importable virus')
        virus_txid = virus_import_parameters[1]
        virus_folder = virus_import_parameters[2]
        generate_fasta(virus_txid, virus_folder, f'{_fasta_target}.fasta')
    elif 'overlap' in action:
        from overlaps import overlaps_controller
        overlaps_controller.run(db_user, db_password, db_port)
    else:
        logger.error(f'the argument {action} is not recognised.\n' + wrong_arguments_message)
except:
    logger.exception('FATAL ERROR')  # this is just to make sure the exception is written to the log file before crashing
    sys.exit(1)
finally:
    if 'import' in action:
        stats_module.check_samples_imported()
    logger.complete()
