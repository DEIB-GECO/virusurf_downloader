import os
from typing import Optional

from Bio import Entrez
from loguru import logger

from data_sources.common_methods_virus import _try_n_times

cached_taxon_id = dict()
cached_taxon_name = dict()

DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_FAILED_PAUSE_SECONDS = 30


def host_taxon_id_from_ncbi_taxon_name(taxon_name: str) -> Optional[int]:
    if not taxon_name:
        return None
    else:
        taxon_id = cached_taxon_id.get(taxon_name.lower())
        if taxon_id == -1:  # -1 means the taxon_id for this taxon name was searched before
            return None
        elif taxon_id is None:
            # fetch it from Entrez API
            def do_api_call():
                global cached_taxon_id
                nonlocal taxon_id
                with Entrez.esearch(db="taxonomy", term=taxon_name, rettype=None,
                                    retmode="xml") as handle:
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
        taxon_name = cached_taxon_name.get(taxon_id)
        if taxon_name == -1:  # -1 means that the taxon_name for this taxon_id was searched before
            return None
        elif taxon_name is None:
            # fetch it from Entrez API
            def do_api_call():
                global cached_taxon_name
                nonlocal taxon_name
                with Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml") as taxon_handle:
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

    """Seems ridiculous, but if you use Efetch on biosample database with the host_sample_id, 
    you download the XML of another host. The only solution I found is to obtain the numerical id
    associated to this host_id, and then use Efetch with that numerical ID.
    See also https://github.com/NCBI-Hackathons/EDirectCookbook/issues/45."""
    def do1():
        with Entrez.esearch(db="biosample", term={host_sample_id}) as handle:
            response = Entrez.read(handle)
        if response['Count'] == '1':
            return int(response['IdList'][0])
        else:
            raise ValueError(f"can't find the biosample numeric id for biosample {host_sample_id}")

    def do2():
        local_file_path = f"{containing_directory}{os.path.sep}{host_sample_id}.xml"
        if not os.path.exists(local_file_path):
            with Entrez.efetch(db="biosample", id=numeric_host_id, retmode="xml") as handle:
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