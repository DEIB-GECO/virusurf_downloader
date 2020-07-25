import os
import time
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Callable

from Bio import Entrez
from lxml import etree

from xml_helper import text_at_node

Entrez.email = "Your.Name.Here@example.org"
from loguru import logger

DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = 30


def _download_virus_taxonomy_as_xml_from_name(containing_directory: str, taxon_name: int) -> str:
    def count_organisms():
        with Entrez.esearch(db="taxonomy", term=f'{taxon_name}', rettype='count', retmode="xml") as handle_1:
            response = Entrez.read(handle_1)
        return int(response['Count'])
    how_many_organisms = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, count_organisms)
    assert how_many_organisms == 1, f'The passed taxon name should match exactly one organism in the NCBI Taxonomy DB but\n' \
                                    f'{taxon_name} matches 0 or more than 1.'

    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{os.path.sep}{taxon_name}.xml"
        if not os.path.exists(destination_file_path):
            # get taxon_id
            with Entrez.esearch(db="taxonomy", term=f'"{taxon_name}"', rettype=None, retmode="xml") as id_search:
                tree: etree.ElementTree = etree.parse(source=id_search, parser=etree.XMLParser(remove_blank_text=True))
                taxon_id = text_at_node(tree, '/eSearchResult/IdList/Id', mandatory=True)

            # download data for taxon_id
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                    open(destination_file_path, 'w') as f:
                for line in handle:
                    f.write(line)
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _download_virus_taxonomy_as_xml(containing_directory: str, taxon_id: int) -> str:
    def count_organisms():
        with Entrez.esearch(db="taxonomy", term=f'{taxon_id}[uid]', rettype='count', retmode="xml") as handle_1:
            response = Entrez.read(handle_1)
        return int(response['Count'])
    how_many_organisms = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, count_organisms)
    assert how_many_organisms == 1, f'The passed taxon id should match exactly one organism in the NCBI Taxonomy DB but\n' \
                                    f'{taxon_id} does not match any.'

    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{os.path.sep}{taxon_id}.xml"
        if not os.path.exists(destination_file_path):
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                    open(destination_file_path, 'w') as f:
                for line in handle:
                    f.write(line)
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)

def _try_n_times(n_times: int, or_wait_secs: int, function: Callable, *args, **kwargs):
    with ThreadPoolExecutor(max_workers=1) as executor:
        try:
            return executor.submit(function, *args, **kwargs).result()
        except IOError as e:
            n_times -= 1
            if n_times > 0:
                logger.error(f'Error while invoking {function.__name__} with args: {args}\nkwargs: {kwargs}. '
                             f'New attempt in {or_wait_secs} secs.'
                             f'Reason of error: {str(type(e))} {e.args}')
                executor.submit(time.sleep, or_wait_secs).result()
                return _try_n_times(n_times, or_wait_secs, function, *args, **kwargs)
            else:
                raise e