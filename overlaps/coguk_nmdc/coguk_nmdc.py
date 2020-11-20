from typing import Optional, List
from overlaps.multi_database_manager import config_db_engine, Session, Sequence, SequencingProject, get_session, rollback, Virus, source_sequences, target_sequences, user_asked_to_commit
from loguru import logger
from tqdm import tqdm
from os.path import sep

# read only
db_name = 'vcm_11_1'
source_database_source = ['COG-UK']
source_name = 'COGUK'
target_database_source = ['NMDC']
target_name = 'NMDC'

total_only_strain_1_to_n = 0
total_strain_plus_length_1_to_n = 0
total_only_strain_1_to_1 = 0
total_strain_plus_length_1_to_1 = 0
output_record = []


def mark_overlaps():
    session = get_session(db_name)
    global total_only_strain_1_to_n, total_strain_plus_length_1_to_n, output_record, total_only_strain_1_to_1, total_strain_plus_length_1_to_1

    try:
        for source_seq in tqdm(source_sequences(session, database_source=source_database_source)):

            only_strain = []
            strain_plus_length = []

            target_seq_query = target_sequences(session,
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
            elif len(only_strain) > 0:
                if len(only_strain) > 1:
                    output_record.append('\t\tWARN\t\t')
                    total_only_strain_1_to_n += 1
                else:
                    total_only_strain_1_to_1 += 1
                output_record.append(f'{source_seq.accession_id} matches with {target_name} strain on {only_strain}')      

        if user_asked_to_commit:
            session.commit()
    except Exception:
        logger.exception("")
        rollback(session)
    finally:
        session.close()

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
    config_db_engine(db_name, db_user, db_password, db_port)
    if not user_asked_to_commit:
        logger.warning('OPERATION WON\'T BE COMMITTED TO THE DB')
    mark_overlaps()

