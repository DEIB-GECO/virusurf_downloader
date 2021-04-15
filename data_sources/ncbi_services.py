from os.path import sep, exists, abspath
from typing import Optional, Tuple, List
from Bio import Entrez
from loguru import logger
from lxml import etree
from xml_helper import text_at_node
from time import sleep
from typing import Callable
from tqdm import tqdm
from http.client import HTTPException
import socket
Entrez.email = "example@mail.com"   # just to silence the warning. Then a correct email can be set later

cached_taxon_id = dict()
cached_taxon_name = dict()
_config_params: Optional[Tuple] = None  # tool name, email, api_key

DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = [10, 30]  # for N attempts, there are N-1 waiting periods
socket.setdefaulttimeout(6)     # sets a timeout on Biopython.Entrez requests (default limit is None = infinite)
                                # affects all socket created from now on


def read_config_params():
    """
    :return: a tuple containing <tool name>, <email>, <api key>
    """
    global _config_params
    if not _config_params:
        config_file_path = f".{sep}e-utils_api_config.csv"
        # create template file if it doesn't exist
        if not exists(config_file_path):
            with open(config_file_path, mode='w') as blank_file:
                blank_file.write("# lines starting with # are comments\n"
                                 "# Write below the following values: tool name,email,api_key")
        # else read it
        values = (None, 'example@email.com', None)  # defaults
        configuration_found = False
        with open(config_file_path, mode='r') as config_file:
            line = config_file.readline()
            while line:
                line = line.strip()
                if not (line.startswith('#') and line != ""):
                    values = line.split(',')
                    assert len(values) == 3, 'E-Utils configuration file contains errors.\n' \
                                             'Expected file contents: tool_name,email,api_key'
                    configuration_found = True
                    break
                line = config_file.readline()
        if not configuration_found:
            logger.warning(f"E-Utils configuration file is missing at path {abspath(config_file_path)}.\n"
                           f"Expected file contents: tool_name,email,api_key\n"
                           f"Press Ctrl+C to correct now or wait 10 seconds for using default parameters (rate-limited):\n"
                           f"{values}")
            try:
                sleep(10)
            except KeyboardInterrupt:
                exit(0)
        Entrez.email = values[1]
        _config_params = values

    return _config_params


def _try_n_times(n_times: int, or_wait_secs: List[int], function: Callable, *args, **kwargs):
    try:
        return function(*args, **kwargs)
    except (RuntimeError, IOError, HTTPException) as e:
        n_times -= 1
        if n_times > 0:
            wait_time = or_wait_secs[-n_times]
            logger.info(f'Error while invoking {function.__name__} with args: {args}\nkwargs: {kwargs}\n'
                        f'\tReason of error: {str(type(e))} {e.args}'
                        f'\tNew attempt in {wait_time}s. Left attempts {n_times}')
            sleep(float(wait_time))
            return _try_n_times(n_times, or_wait_secs, function, *args, **kwargs)
        else:
            raise e
            

def host_taxon_id_from_ncbi_taxon_name(taxon_name: str) -> Optional[int]:
    if not taxon_name:
        return None
    else:
        global cached_taxon_id
        taxon_id = cached_taxon_id.get(taxon_name.lower())
        if taxon_id == -1:  # -1 means the taxon_id for this taxon name was searched before
            return None
        elif taxon_id is None:
            # fetch it from Entrez API
            entrez_config = read_config_params()
            def do_api_call():
                nonlocal taxon_id
                with Entrez.esearch(db="taxonomy", term=taxon_name, rettype=None, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle:
                    response = Entrez.read(handle)
                    if response['Count'] == '1':
                        taxon_id = int(response['IdList'][0])
                        cached_taxon_id[taxon_name] = taxon_id
                    else:
                        logger.warning(f'can\'t find the taxon id for taxon name {taxon_name}')
                        cached_taxon_id[
                            taxon_name] = -1  # save -1 in cache to distinguish from non cached taxon_ids

            try:
                _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do_api_call)
            except:
                logger.exception(f'Impossible to fetch the taxon_id corresponding to {taxon_name}.')

        return taxon_id


def host_taxon_name_from_ncbi_taxon_id(taxon_id: int) -> Optional[str]:
    if not taxon_id:
        return None
    else:
        global cached_taxon_name
        taxon_name = cached_taxon_name.get(taxon_id)
        if taxon_name == -1:  # -1 means that the taxon_name for this taxon_id was searched before
            return None
        elif taxon_name is None:
            entrez_config = read_config_params()
            # fetch it from Entrez API
            def do_api_call():
                nonlocal taxon_name
                with Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as taxon_handle:
                    response = Entrez.read(taxon_handle)  # response is an array of taxons
                    if len(response) > 0:
                        taxon_name = str(response[0]['ScientificName'])  # parse from Bio.Entrez.Parser.StringElement
                        cached_taxon_name[taxon_id] = taxon_name.lower()
                    else:
                        logger.warning(f'can\'t find the taxon_name for taxon_id {taxon_id}')
                        cached_taxon_name[taxon_id] = -1  # save -1 in cache to distinguish from non cached taxon_ids

            try:
                _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do_api_call)
            except:
                logger.exception(f'Exception occurred while fetching the taxon_name corresponding to id {taxon_id}.')

        return taxon_name


def download_or_get_ncbi_host_sample_as_xml(containing_directory: str, host_sample_id: str) -> str:
    """
    :param containing_directory: directory where the file will be downloaded and cached
    :param host_sample_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    entrez_config = read_config_params()

    """Seems ridiculous, but if you use Efetch on biosample database with the host_sample_id, 
    you download the XML of another host. The only solution I found is to obtain the numerical id
    associated to this host_id, and then use Efetch with that numerical ID.
    See also https://github.com/NCBI-Hackathons/EDirectCookbook/issues/45."""
    def do1():
        with Entrez.esearch(db="biosample", term={host_sample_id}, tool=entrez_config[0], email=entrez_config[1],
                            api_key=entrez_config[2]) as handle:
            response = Entrez.read(handle)
        if response['Count'] == '1':
            return int(response['IdList'][0])
        else:
            raise ValueError(f"can't find the biosample numeric id for biosample {host_sample_id}")

    def do2():
        local_file_path = f"{containing_directory}{sep}{host_sample_id}.xml"
        if not exists(local_file_path):
            with Entrez.efetch(db="biosample", id=numeric_host_id, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle:
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
    numeric_host_id = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do1)
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do2)


def download_ncbi_taxonomy_as_xml_from_name(containing_directory: str, taxon_name: str) -> str:
    entrez_config = read_config_params()
    def count_organisms():
        with Entrez.esearch(db="taxonomy", term=f'{taxon_name}', rettype='count', retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle_1:
            response = Entrez.read(handle_1)
        return int(response['Count'])
    how_many_organisms = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, count_organisms)
    assert how_many_organisms == 1, f'The passed taxon name should match exactly one organism in the NCBI Taxonomy DB but\n' \
                                    f'{taxon_name} matches 0 or more than 1.'

    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{sep}{taxon_name}.xml"
        if not exists(destination_file_path):
            # get taxon_id
            with Entrez.esearch(db="taxonomy", term=f'"{taxon_name}"', rettype=None, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as id_search:
                tree: etree.ElementTree = etree.parse(source=id_search, parser=etree.XMLParser(remove_blank_text=True))
                taxon_id = text_at_node(tree, '/eSearchResult/IdList/Id', mandatory=True)

            # download data for taxon_id
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle:
                with open(destination_file_path, 'w') as f:
                    f.write(handle.read())
        return destination_file_path
    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)


def download_ncbi_taxonomy_as_xml(containing_directory: str, taxon_id: int) -> str:
    entrez_config = read_config_params()
    def count_organisms():
        with Entrez.esearch(db="taxonomy", term=f'{taxon_id}[uid]', rettype='count', retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle_1:
            response = Entrez.read(handle_1)
        return int(response['Count'])
    how_many_organisms = _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, count_organisms)
    assert how_many_organisms == 1, f'The passed taxon id should match exactly one organism in the NCBI Taxonomy DB but\n' \
                                    f'{taxon_id} does not match any.'

    def do():
        # download and write the taxonomy tree for the taxon id
        destination_file_path = f"{containing_directory}{sep}{taxon_id}.xml"
        if not exists(destination_file_path):
            with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle:
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
        local_file_path = f"{containing_directory}{sep}{sample_accession_id}.xml"
        if not exists(local_file_path):
            entrez_config = read_config_params()
            with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml", tool=entrez_config[0],
                email=entrez_config[1], api_key=entrez_config[2]) as handle:
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


def get_samples_accession_ids(samples_query: str) -> List[int]:
    logger.info(f'getting accession ids of samples...')

    def do():
        entrez_config = read_config_params()
        # DO PAGINATION
        # total number of sequences
        with Entrez.esearch(db="nuccore", term=f"{samples_query}", rettype='count', tool=entrez_config[0],
                            email=entrez_config[1], api_key=entrez_config[2]) as handle:
            response = Entrez.read(handle)
            total_records = int(response['Count'])
        # get pages
        accessions_ids = list()
        # noinspection PyPep8Naming
        RECORDS_PER_PAGE = 5000
        page_number = 0
        import time
        with tqdm(total=total_records) as progress_bar:
            while total_records > page_number * RECORDS_PER_PAGE:
                with Entrez.esearch(db="nuccore",
                                    term=f"{samples_query}",
                                    retmax=RECORDS_PER_PAGE, retstart=page_number * RECORDS_PER_PAGE) as handle:
                    response = Entrez.read(handle)
                for x in response['IdList']:
                    accessions_ids.append(int(x))
                page_number += 1
                time.sleep(2)
                progress_bar.update(
                    RECORDS_PER_PAGE if page_number * RECORDS_PER_PAGE < total_records else total_records - (
                                (page_number - 1) * RECORDS_PER_PAGE))
        if len(accessions_ids) != total_records:
            raise IOError('Some of the accession ids were not correctly downloaded')
        return accessions_ids

    return _try_n_times(DOWNLOAD_ATTEMPTS, DOWNLOAD_FAILED_PAUSE_SECONDS, do)