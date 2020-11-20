"""
Created by tomalf2 on nov, 2020.
"""
from os.path import sep, abspath
from loguru import logger

file_path = f'.{sep}overlaps{sep}genbank_gisaid{sep}genbank_gisaid_overlaps.txt'

logger.info(f'Reading file: {abspath(file_path)}')

def read_source_matching_ids():
    ids = []
    with open(file_path, mode='r') as f:
        for line in f.readlines():
            if line.startswith('\t'):
                continue
            elif line.startswith('TOTAL'):
                break
            else:
                end_idx = line.index(' ')
                ids.append(line[:end_idx])

    print(ids)


read_source_matching_ids()