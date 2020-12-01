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


def download_ncbi_taxonomy_as_xml_from_name(containing_directory: str, taxon_name: str) -> str:
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
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle:
                with open(destination_file_path, 'w') as f:
                    f.write(handle.read())
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def download_ncbi_taxonomy_as_xml(containing_directory: str, taxon_id: int) -> str:
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
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle:
                try:
                    with open(destination_file_path, 'w') as f:
                        f.write(handle.read())
                except TypeError:
                    with open(destination_file_path, 'wb') as f:
                        f.write(handle.read())
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def download_or_get_ncbi_sample_as_xml(containing_directory: str, sample_accession_id: int) -> str:
    """
    :param containing_directory: directory where the file will be downloaded and cached
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    def do():
        local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
        if not os.path.exists(local_file_path):
            with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml") as handle:
                a = handle.read()
                try:
                    with open(local_file_path, 'w') as f:
                        f.write(a)
                except TypeError:
                    # sometimes the source handle returns bytes instead of chars, so let's try writing bytes
                    try:
                        with open(local_file_path, 'wb') as f:
                            f.write(a)
                    except Exception as e:
                        if a:
                            logger.error(f'Content of EntrezAPI was:\n{a[:50]}... (only 1st 50 chars displayed')
                        else:
                            logger.error(f'Content of EntrezAPI is probably empty.')
                        raise e
        return local_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def _try_n_times(n_times: int, or_wait_secs: int, function: Callable, *args, **kwargs):
    with ThreadPoolExecutor(max_workers=1) as executor:
        try:
            return executor.submit(function, *args, **kwargs).result()
        except (RuntimeError, IOError) as e:
            n_times -= 1
            if n_times > 0:
                logger.info(f'Error while invoking {function.__name__} with args: {args}\nkwargs: {kwargs}\n'
                             f'\tReason of error: {str(type(e))} {e.args}'
                             f'\tNew attempt in {or_wait_secs}s. Left attempts {n_times}')
                executor.submit(time.sleep, or_wait_secs).result()
                return _try_n_times(n_times, or_wait_secs, function, *args, **kwargs)
            else:
                raise e