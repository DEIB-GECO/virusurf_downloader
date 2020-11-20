from typing import Optional, List
from overlaps.multi_database_manager import config_db_engine, Session, Sequence, SequencingProject, get_session, rollback, Virus, source_sequences, target_sequences, user_asked_to_commit
from sqlalchemy import func, or_
from loguru import logger
from tqdm import tqdm
from os.path import sep

# read only values
source_db_name = 'vcm_gisaid_12'
source_name = 'GISAID'
source_database_source = None
target_db_name = 'vcm_11_1'
target_name = 'NMDC'
target_database_source = ['NMDC']


total_only_strain_1_to_n = 0
total_strain_plus_length_1_to_n = 0
total_only_strain_1_to_1 = 0
total_strain_plus_length_1_to_1 = 0
output_record = []


def mark_overlaps():
    source_session = get_session(source_db_name)
    target_session = get_session(target_db_name)
    global total_only_strain_1_to_n, total_strain_plus_length_1_to_n, output_record, total_only_strain_1_to_1, total_strain_plus_length_1_to_1

    try:
        for source_seq in tqdm(source_sequences(source_session)):

            only_strain = []
            strain_plus_length = []

            target_seq_query = target_sequences(target_session,
                                                matching_strain=source_seq.strain_name,
                                                database_source=target_database_source)

            for target_seq in target_seq_query:
                if target_seq.length == source_seq.length:
                    strain_plus_length.append(target_seq.accession_id)
                else:
                    only_strain.append(target_seq.accession_id)

            if len(strain_plus_length) > 0:
                if len(strain_plus_length) > 1:
                    output_record.append('\t\tWARN\t\t')
                    total_strain_plus_length_1_to_n += 1
                else:
                    total_strain_plus_length_1_to_1 += 1
                output_record.append(f'{source_seq.accession_id} matches with {target_name} strain+length on {strain_plus_length}')
                source_seq.gisaid_only = False
            elif len(only_strain) > 0:
                if len(only_strain) > 1:
                    output_record.append('\t\tWARN\t\t')
                    total_only_strain_1_to_n += 1
                else:
                    total_only_strain_1_to_1 += 1
                output_record.append(f'{source_seq.accession_id} matches with {target_name} strain on {only_strain}')
                source_seq.gisaid_only = False

        if user_asked_to_commit:
            source_session.commit()
    except Exception:
        logger.exception("")
        rollback(source_session)
        rollback(target_session)
    finally:
        source_session.close()
        target_session.close()

        totals_string = f'TOTALS:\n' \
                        f'1-1 MATCHES: strain+length: {total_strain_plus_length_1_to_1} -- only strain {total_only_strain_1_to_1}\n' \
                        f'1-N MATCHES: strain+length: {total_strain_plus_length_1_to_n} -- only strain {total_only_strain_1_to_n} (search "WARN" to find \'em)\n'
        logger.info(totals_string)

        output_path = f'.{sep}overlaps{sep}{source_name}_{target_name}{sep}{source_name}_{target_name}_overlaps.txt'.lower()
        with open(file=output_path, mode='w') as w:
            for line in output_record:
                w.write(line + "\n")
            w.write(totals_string)


def run(db_user, db_password, db_port):
    config_db_engine(source_db_name, db_user, db_password, db_port)
    config_db_engine(target_db_name, db_user, db_password, db_port)
    if not user_asked_to_commit:
        logger.warning('OPERATION WON\'T BE COMMITTED TO THE DB')
    mark_overlaps()

