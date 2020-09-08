import sys
from loguru import logger
import database_tom
from data_sources.ncbi_any_virus.ncbi_importer import prepared_parameters, import_samples_into_vcm
import data_sources.nmdc.procedure as nmdc
from Bio import Entrez
from tqdm import tqdm
import stats_module
import warnings

Entrez.email = "Your.Name.Here@example.org"


#   #################################       PROGRAM ARGUMENTS   ##########################
wrong_arguments_message = 'The module main.py expects the following arguments:' \
                          'db_name, recreate_db?, db_user, db_password, db_port, source_to_import\n' \
                          'source_to_import accept values:\n' \
                          'cog-uk\n' \
                          'gisaid\n'
for p in prepared_parameters.keys():
    wrong_arguments_message += f'{p}\n'
wrong_arguments_message += 'just_make_indexes\n'
wrong_arguments_message += 'create_views\n'
wrong_arguments_message += 'disambiguate_chimera_sequences\n'
wrong_arguments_message += 'epitopes\n'
wrong_arguments_message += 'If you want to import virus samples, you can optionally specify a range of samples to import as <min> (included) <max> (excluded).\n'
wrong_arguments_message += 'If you want to import epitopes, you must specify a source name between the above ones.\n'
# noinspection PyBroadException
try:
    db_name = sys.argv[1]
    _db_recreate = sys.argv[2].lower()
    db_recreate = True if _db_recreate == 'yes' or _db_recreate == 'true' else False
    db_user = sys.argv[3]
    db_password = sys.argv[4]
    db_port = sys.argv[5]
    source = sys.argv[6].lower()
    try:
        _from = int(sys.argv[7])
        to = int(sys.argv[8])
        if not _from < to:
            logger.error('Optional parameter <min> must be less than <max>')
            sys.exit(1)
    except Exception:
        _from = None
        to = None
    if 'epitope' in source:
        _epitope_target = sys.argv[7]
        if not _epitope_target:
            logger.error('No epitope target specified. Type the target virus after "epitopes"')
            raise Exception()
except Exception:
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
logger.remove()  # removes default logger to stderr with level DEBUG
# on console print from level INFO on
logger.add(sink=lambda msg: tqdm.write(msg, end=''),
           level='TRACE',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=True,
           backtrace=True,
           diagnose=True,
           enqueue=True)
# log to file any message of any security level
logger.add("./logs/log_"+source+"_{time}.log",
           level='TRACE',
           rotation='100 MB',
           format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
                  "<level>{level: <8}</level> | "
                  "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
           colorize=False,
           backtrace=True,
           diagnose=True,
           enqueue=True)
# redirect warnings
def customwarn(message, category, filename, lineno, file=None, line=None):
    logger.warning(warnings.formatwarning(message, category, filename, lineno))
warnings.showwarning = customwarn

logger.info(f"main.py {' '.join(sys.argv[1:])}")

#   ###################################     FILL DB WITH VIRUS SEQUENCES    ###############
# init database
database_tom.config_db_engine(db_name, db_user, db_password, db_port, recreate_db_from_scratch=db_recreate)

#   ###################################     MAIN OPERATION       ###############
try:
    if 'index' in source:
        database_tom.create_indexes()
    elif 'view' in source:
        database_tom.create_views()
    elif 'chimera_sequence' in source:
        database_tom.disambiguate_chimera_sequences()
    elif 'epitope' in source:
        from epitopes import import_epitopes
        import_epitopes()
    elif source in ['coguk', 'cog-uk']:
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
except:
    logger.exception('FATAL ERROR')  # this is just to make sure the exception is written to the log file before crashing
finally:
    stats_module.check_samples_imported()
    logger.complete()
