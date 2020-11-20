"""
Created by tomalf2 on nov, 2020.
"""
from os.path import abspath
from database_tom import try_py_function, Sequence
from sqlalchemy import func
from loguru import logger
from typing import List
from tqdm import tqdm


def load(lineages_tsv_file_path):

    def check_lineage_are_assigned_in_db(session):
        accession_ids = [x[0] for x in read_lineages()]
        not_null_lineages = session\
            .query(func.count(Sequence.sequence_id))\
            .filter(Sequence.accession_id.in_(accession_ids),
                    Sequence.lineage.isnot(None)).scalar()
        logger.info(f'Number of sequences with not-null lineages {not_null_lineages}')

    def update_with_lineages(session):
        for accession_id, lineage in tqdm(read_lineages()):
            sequence = session\
                .query(Sequence)\
                .filter(Sequence.accession_id == accession_id).one_or_none()
            sequence.lineage = lineage

    def read_lineages():
        with open(lineages_tsv_file_path, mode='r') as lineages_file:
            lineages_file.readline()    # skip header line
            for line in lineages_file:
                # extract interesting values
                first_comma_idx = line.find(',')
                accession_id = line[:first_comma_idx]
                second_comma_idx = line.find(',', first_comma_idx + 1)
                lineage = line[first_comma_idx + 1:second_comma_idx]
                if lineage == 'None':
                    continue
                yield accession_id, lineage

    logger.info('BEFORE UPDATING LINEAGES')
    try_py_function(
        check_lineage_are_assigned_in_db
    )

    logger.info('UPDATING LINEAGES...')
    try_py_function(
        update_with_lineages
    )

    logger.info('AFTER UPDATING LINEAGES')
    try_py_function(
        check_lineage_are_assigned_in_db
    )