"""
Created by tomalf2 on ott, 2020.
"""
from typing import Generator, Tuple
from sqlalchemy import func
from tqdm import tqdm
from db_config import read_db_import_configuration as import_config, database
from locations import get_local_folder_for, FileType
from epitopes import virus_database_id
from os.path import abspath
from loguru import logger


def generate_fasta(virus_taxon_id: int, virus_folder_name: str, generated_file_name: str,
                   only_null_lineages: bool = False) -> str:
    """
    Generates a multi fasta file containing all the sequences of the given virus.

    :return: the absolute path to the generated fasta file.
    """
    db_params: dict = import_config.get_database_config_params()
    database.config_db_engine(db_params["db_name"], db_params["db_user"], db_params["db_psw"], db_params["db_port"])

    virus_db_id = virus_database_id(virus_taxon_id)

    if virus_db_id is None:
        raise Exception('Before running this algorithm, create the '
                        f'virus associated with taxon {virus_taxon_id}')

    def get_acc_ids_and_sequences_from_db(session: database.Session) -> Generator[Tuple, None, None]:
        query = session.query(database.Sequence.accession_id,
                                     database.NucleotideSequence.nucleotide_sequence)\
            .filter(database.Sequence.virus_id == virus_db_id,
                    database.Sequence.sequence_id == database.NucleotideSequence.sequence_id)
        if only_null_lineages:
            query = query.filter(database.Sequence.lineage == None)
        for pair in query.all():
            yield pair[0], pair[1]

    def get_total_acc_ids_from_db(session: database.Session) -> int:
        query = session.query(func.count(database.Sequence.accession_id)) \
            .filter(database.Sequence.virus_id == virus_db_id)
        if only_null_lineages:
            query = query.filter(database.Sequence.lineage == None)
        return query.first()[0]

    target_file_path = get_local_folder_for(virus_folder_name, FileType.Fasta) + generated_file_name
    logger.info(f"Generating fasta...")
    with open(file=target_file_path, mode='w') as file:
        total_count = database.try_py_function(get_total_acc_ids_from_db)
        print(total_count)
        data = database.try_py_function(get_acc_ids_and_sequences_from_db)
        progress = tqdm(total=total_count)
        if total_count > 0:
            first = next(data)
            file.write(f'>{first[0]}\n')
            file.write(first[1])
            progress.update()
        for acc_id_and_sequences in data:
            file.write(f'\n>{acc_id_and_sequences[0]}\n')
            file.write(acc_id_and_sequences[1])
            progress.update()
    target_file_path = abspath(target_file_path)
    logger.info(f"Fasta file generated at {target_file_path}")
    return target_file_path
