"""
Created by tomalf2 on ott, 2020.
"""
from typing import Generator, Tuple, Collection
from sqlalchemy import func
from tqdm import tqdm
from db_config import read_db_import_configuration as import_config, database
from locations import get_local_folder_for, FileType
from epitopes import virus_database_id
from os.path import abspath
from loguru import logger
from db_config.database import NucleotideSequence, SequencingProject, Sequence, Session


def generate_fasta(virus_taxon_id: int, virus_folder_name: str, generated_file_name: str,
                   only_null_lineages: bool = False, data_source: Collection[str] = None) -> str:
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

    def get_acc_ids_and_sequences_from_db(session: Session) -> Generator[Tuple, None, None]:
        tables_in_from = [Sequence, NucleotideSequence]
        if data_source:
            tables_in_from.append(SequencingProject)

        query = session.query(Sequence.accession_id,
                              NucleotideSequence.nucleotide_sequence)\
            .select_from(*tables_in_from)\
            .filter(Sequence.virus_id == virus_db_id,
                    Sequence.sequence_id == NucleotideSequence.sequence_id)
        if only_null_lineages:
            query = query.filter(Sequence.lineage == None)
        if data_source:
            query = query.filter(
                Sequence.sequencing_project_id == SequencingProject.sequencing_project_id,
                SequencingProject.database_source.in_(data_source)
            )
        for pair in query.all():
            yield pair[0], pair[1]

    def get_total_acc_ids_from_db(session: Session) -> int:
        tables_in_from = [Sequence]
        if data_source:
            tables_in_from.append(SequencingProject)

        query = session.query(func.count(Sequence.accession_id))\
            .select_from(*tables_in_from)\
            .filter(Sequence.virus_id == virus_db_id)
        if only_null_lineages:
            query = query.filter(Sequence.lineage == None)
        if data_source:
            query = query.filter(
                Sequence.sequencing_project_id == SequencingProject.sequencing_project_id,
                SequencingProject.database_source.in_(data_source)
            )
        return query.first()[0]

    target_file_path = get_local_folder_for(virus_folder_name, FileType.Fasta) + generated_file_name
    logger.info(f"Generating fasta...")
    with open(file=target_file_path, mode='w') as file:
        total_count = database.try_py_function(get_total_acc_ids_from_db)
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


if __name__ == '__main__':
    import sys
    from time import sleep

    print('Arguments:\n'
          '- db_name\n'
          '- virus_taxon_id\n'
          '- virus dir name\n'
          '- output fasta name\n'
          '- inlcude only sequences with null lineage ?\n'
          '- all_except_coguk ? ')

    db_name = sys.argv[1]
    v_taxon_id = int(sys.argv[2])
    v_dir_name = str(sys.argv[3])
    generated_file_name = str(sys.argv[4])
    include_only_null_lineages = True if str(sys.argv[5]).lower() in ('true', 'yes', 'only_new') else False
    data_sources = ('GenBank', 'NMDC', 'RefSeq') if sys.argv[6].lower() in ('true', 'yes') else None

    print('ARGUMENTS:')
    print(db_name, v_taxon_id, v_dir_name, generated_file_name, include_only_null_lineages, data_sources)
    sleep(10)

    import_config.set_db_name(db_name)
    generate_fasta(v_taxon_id, v_dir_name, generated_file_name, include_only_null_lineages, data_sources)
