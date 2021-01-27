from typing import Optional, List
from overlaps.multi_database_manager import config_db_engine, Session, Sequence, SequencingProject, get_session, \
    rollback, Virus, source_sequences, target_sequences, user_asked_to_commit, insert_overlaps_in_db, \
    cleanup_overlap_tables
from loguru import logger
from tqdm import tqdm
from os.path import sep
from datetime import date
import db_config.read_db_overlaps_configuration as db_config

# read only
db_name = ''
source_database_source = ['GenBank', 'RefSeq']
source_name = 'GenBank'
target_database_source = ['NMDC']
target_name = 'NMDC'

total_only_strain_1_to_n = 0
total_strain_plus_length_1_to_n = 0
total_only_strain_1_to_1 = 0
total_strain_plus_length_1_to_1 = 0
output_record = []


def mark_overlaps():
    session = get_session(db_name)
    cleanup_overlap_tables(session)
    global total_only_strain_1_to_n, total_strain_plus_length_1_to_n, output_record, total_only_strain_1_to_1, total_strain_plus_length_1_to_1

    try:
        count_source_seq = source_sequences(session=session,
                                            database_source=source_database_source,
                                            virus_taxon_name='Severe acute respiratory syndrome coronavirus 2',
                                            for_overlaps_with_target_source=target_name,
                                            count_only=True)
        for source_seq in tqdm(total=count_source_seq,
                               iterable=source_sequences(session=session,
                                                         database_source=source_database_source,
                                                         virus_taxon_name='Severe acute respiratory syndrome coronavirus 2',
                                                         for_overlaps_with_target_source=target_name)):

            only_strain = []
            strain_plus_length = []

            target_seq_query = target_sequences(session=session,
                                                matching_strain=source_seq.strain_name,
                                                database_source=target_database_source)

            for target_seq in target_seq_query:
                if target_seq.length == source_seq.length:
                    strain_plus_length.append(target_seq)
                else:
                    only_strain.append(target_seq)

            if len(strain_plus_length) > 0:
                if len(strain_plus_length) > 1:
                    output_record.append('\t\tWARN\t\t')
                    total_strain_plus_length_1_to_n += 1
                else:
                    total_strain_plus_length_1_to_1 += 1
                acc_ids = [s.accession_id for s in strain_plus_length]
                output_record.append(f'{source_seq.accession_id} matches with {target_name} strain+length on {acc_ids}')
                insert_overlaps_in_db(session, session, source_seq, strain_plus_length, source_name, target_name)
            elif len(only_strain) > 0:
                if len(only_strain) > 1:
                    output_record.append('\t\tWARN\t\t')
                    total_only_strain_1_to_n += 1
                else:
                    total_only_strain_1_to_1 += 1
                acc_ids = [s.accession_id for s in only_strain]
                output_record.append(f'{source_seq.accession_id} matches with {target_name} strain on {acc_ids}')
                insert_overlaps_in_db(session, session, source_seq, only_strain, source_name, target_name)

        if user_asked_to_commit:
            session.commit()
    except KeyboardInterrupt:
        rollback(session)
        output_record.append("COMPUTATION INTERRUPTED. TOTALS MAY BE INCOMPLETE !!")
    except Exception as e:
        logger.exception("")
        rollback(session)
        output_record.append("COMPUTATION INTERRUPTED. TOTALS MAY BE INCOMPLETE !!")
    finally:
        session.close()

        totals_string = f'TOTALS:\n' \
                        f'1-1 MATCHES: strain+length: {total_strain_plus_length_1_to_1} -- only strain {total_only_strain_1_to_1}\n' \
                        f'1-N MATCHES: strain+length: {total_strain_plus_length_1_to_n} -- only strain {total_only_strain_1_to_n} (search "WARN" to find \'em)\n'
        logger.info(totals_string)

        output_path = f'.{sep}overlaps{sep}{source_name}_{target_name}{sep}'
        output_path += f'{source_name}@{db_name}_overlapping_{target_name}@{db_name}'
        output_path += f'_{date.today().strftime("%Y-%b-%d")}.txt'
        output_path = output_path.lower()
        with open(file=output_path, mode='w') as w:
            for line in output_record:
                w.write(line + "\n")
            w.write(totals_string)


def run():
    db = db_config.get_import_params_for("genbank")
    global db_name
    db_name = db["db_name"]
    config_db_engine(db["db_name"], db["db_user"], db["db_psw"], db["db_port"])
    if not user_asked_to_commit:
        logger.warning('OPERATION WON\'T BE COMMITTED TO THE DB')
    mark_overlaps()


