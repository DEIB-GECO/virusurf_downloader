import sys
from loguru import logger
from data_sources.ncbi_any_virus.ncbi_importer import import_samples_into_vcm
from data_sources.ncbi_any_virus.settings import known_settings
import data_sources.nmdc.procedure as nmdc
import stats_module
from logger_settings import setup_logger
from db_config import read_db_import_configuration as import_config, database_tom

log_file_keyword = ""

#   #################################       PROGRAM ARGUMENTS   ##########################
ncbi_virus_names = sorted(known_settings.keys())
wrong_arguments_message = "Accepted parameters:\n" \
                          "\timport <db_name> coguk|gisaid|" + "|".join(ncbi_virus_names) + "\n" \
                          "\tepitopes <db_name> " + "|".join(ncbi_virus_names) + "\n" \
                          "\toverlaps <source1>_<source2> <commit to the DB ?>\n" \
                          "\tlineages <db_name> <lineages TSV file path>"
wrong_arguments_message += 'When action is "import", the source name can be optionally followed by a range of samples to import as <min> (included) <max> (excluded).\n'

try:
    action = sys.argv[1].lower()
    if action in ['import', 'epitope', 'lineage', 'fasta', 'sql']:     # get database name if action is not 'overlaps'
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
    elif 'epitope' in action:
        _epitope_target = sys.argv[3]
        log_file_keyword = f"epi_{_epitope_target}"
    elif 'lineage' in action:
        _lineages_tsv_file_path = sys.argv[3]
        log_file_keyword = f"lineages"
    elif 'fasta' in action:
        _fasta_target = sys.argv[3]
        log_file_keyword = f"fasta_{_fasta_target}"
    elif 'overlap' in action:
        _overlap_target = sys.argv[2]
        log_file_keyword = f"overlaps_{_overlap_target}"
    else:
        log_file_keyword = action
except IndexError:
    logger.error(wrong_arguments_message)
    sys.exit(1)

#   ###################################      SETUP LOGGER    ##############################
setup_logger(log_file_keyword)


#   ###################################     PERFORM <action>       ###############
try:
    if 'index' in action:
        database_tom.create_indexes()
    elif 'view' in action:
        database_tom.create_views()
    elif 'chimera_sequence' in action:
        database_tom.disambiguate_chimera_sequences()
    elif 'epitope' in action:
        from epitopes import import_epitopes
        # noinspection PyUnboundLocalVariable
        virus_import_parameters = known_settings.get(_epitope_target)
        if not virus_import_parameters:
            raise ValueError(f'{_epitope_target} is not recognised as an importable virus')
        virus_txid = virus_import_parameters[1]
        import_epitopes(virus_txid)
    elif 'lineage' in action:
        from load_lineages import load
        # noinspection PyUnboundLocalVariable
        load(_lineages_tsv_file_path)
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
    elif 'fasta' in action:
        from generate_fasta import generate_fasta
        # noinspection PyUnboundLocalVariable
        virus_import_parameters = known_settings.get(_fasta_target)
        if not virus_import_parameters:
            raise ValueError(f'{_fasta_target} is not recognised as an importable virus')
        virus_txid = virus_import_parameters["virus_taxon_id"]
        virus_folder = virus_import_parameters["generated_dir_name"]
        generate_fasta(virus_txid, virus_folder, f'{_fasta_target}.fasta')
    elif 'overlap' in action:
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
