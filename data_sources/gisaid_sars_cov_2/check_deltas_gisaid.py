from overlaps.multi_database_manager import config_db_engine, source_sequences, target_sequences, get_session, Sequence
from sqlalchemy import func
from os.path import sep
from tqdm import tqdm
from locations import get_local_folder_for, FileType

# read only values
source_db_name = 'vcm_gisaid_11'
source_name = 'GISAID_11'
target_db_name = 'vcm_gisaid_12'
target_name = 'GISAID_12'

# source_seq_file_path = f'.{sep}overlaps{sep}check_deltas_gisaid{sep}vcm_gisaid_11_seq.csv'
# target_seq_file_path = f'.{sep}overlaps{sep}check_deltas_gisaid{sep}vcm_gisaid_12_seq.csv'
source_seq_file_path = f"{get_local_folder_for('GISAID', FileType.Other)}vcm_gisaid_11_seq.csv"
source_seq_file_path = f"{get_local_folder_for('GISAID', FileType.Other)}vcm_gisaid_12_seq.csv"


def write_comparable_infos_to_csv(db_name, file_path):
    session = get_session(db_name)
    total_rows = None
    rows_written = 0
    try:
        total_rows = session.query(func.count(Sequence.sequence_id)).scalar()

        with open(file_path, mode='w') as file:
            for source_seq in tqdm(session.query(Sequence) \
                .filter(Sequence.strain_name.isnot(None)) \
                .order_by(Sequence.accession_id).yield_per(100)):

                # interesting atttributes
                string = ','.join([
                    source_seq.accession_id,
                    source_seq.strain_name,
                    str(source_seq.length)
                ])

                file.write(string+'\n')
                rows_written += 1
    except Exception as e:
        print(e)
    finally:
        session.close()

    print(f'written {rows_written} / {total_rows}. Missing rows had no strain_name.')



def compare_infos(smaller_file_path, larger_file_path):
    with open(smaller_file_path, mode='r') as smaller_f, open(larger_file_path, mode='r') as larger_f:
        last_line_read_from_larger_file = 0

        matching_ids_found = 0
        matching_everything = 0
        matching_length_only = 0
        matching_strain_only = 0
        total_sequences_searched = 0
        total_sequences_removed = 0

        try:
            searched_seq = smaller_f.readline()
            matching_seq = larger_f.readline()
            while searched_seq:
                searched_seq = searched_seq.rstrip()

                ac_id, strain, length = searched_seq.split(',')
                matching_id = False
                while matching_seq and not matching_id:
                    last_line_read_from_larger_file += 1
                    matching_seq = matching_seq.rstrip()

                    mat_ac_id, mat_strain, mat_length = matching_seq.split(',')
                    # check to not having surpassed the source acc_id
                    if mat_ac_id > ac_id:
                        # source id is missing in target file
                        total_sequences_removed += 1
                        break   # skip update instruction otherwise all ids from now on will be grater than source ids
                    elif ac_id == mat_ac_id:
                        matching_id = True
                        matching_ids_found += 1
                        if strain == mat_strain and length == mat_length:
                            matching_everything += 1
                        elif length == mat_length:
                            matching_length_only += 1
                        elif strain == mat_strain:
                            matching_strain_only += 1
                    matching_seq = larger_f.readline()

                searched_seq = smaller_f.readline()
                total_sequences_searched += 1
        except Exception as e:
            print(e)
        finally:
            print(
                f'sequences evaluated from smaller file: {total_sequences_searched}\n'
                f'matching accession_ids found: {matching_ids_found}\n'
                f'\twith strain+length equal: {matching_everything}\n'
                f'\twith only length equal: {matching_length_only}\n'
                f'\twith only strain equal: {matching_strain_only}\n'
                f'removed sequences: {total_sequences_removed}\n'
            )
            print(f'lines read from {target_name}: {last_line_read_from_larger_file}')



def find_sequences_removed_from_remote():
    session_11 = get_session(source_db_name)
    old_accession_ids = set(session_11.query(Sequence.accession_id).all())
    session_11.close()
    session_12 = get_session(target_db_name)
    new_accession_ids = set(session_12.query(Sequence.accession_id).all())
    session_12.close()
    removed_ids = old_accession_ids - new_accession_ids
    print(f'Number of sequences present in vcm_gisaid_11 and missing in vcm_gisaid_12: {len(removed_ids)}')



def run(db_user, db_password, db_port):

    config_db_engine(source_db_name, db_user, db_password, db_port)
    # write_comparable_infos_to_csv(source_db_name, source_seq_file_path)
    config_db_engine(target_db_name, db_user, db_password, 5432)
    # write_comparable_infos_to_csv(target_db_name, target_seq_file_path)

    # compare_infos(source_seq_file_path, target_seq_file_path)
    find_sequences_removed_from_remote()

# sequences evaluated from smaller file: 109922
# matching accession_ids found: 109871
#         with strain+length equal: 109221
#         with only length equal: 633
#         with only strain equal: 11
# removed sequences: 51
#
# lines read from GISAID_12: 109958

