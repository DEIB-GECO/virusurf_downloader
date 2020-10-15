from typing import Optional

from Bio import Entrez
from loguru import logger

cached_taxon_id = {}
cached_taxon_name = {}


def host_taxon_id_from_ncbi_taxon_name(taxon_name: str) -> Optional[int]:
    if not taxon_name:
        return None
    else:
        global cached_taxon_id
        taxon_id = cached_taxon_id.get(taxon_name.lower())
        if taxon_id == -1:  # -1 means the taxon_id for this taxon name was searched before
            return None
        elif taxon_id is None:
            try:
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
                        taxon_id = None
            except:
                logger.exception(f'Exception occurred while fetching the taxon_id corresponding to {taxon_name}.')
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
            try:
                with Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml") as taxon_handle:
                    response = Entrez.read(taxon_handle)  # response is an array of taxons
                    if len(response) > 0:
                        taxon_name = str(response[0]['ScientificName'])  # parse from Bio.Entrez.Parser.StringElement
                        cached_taxon_name[taxon_id] = taxon_name.lower()
                    else:
                        logger.warning(f'can\'t find the taxon_name for taxon_id {taxon_id}')
                        cached_taxon_name[taxon_id] = -1  # save -1 in cache to distinguish from non cached taxon_ids
                        taxon_name = None
            except:
                logger.exception(f'Exception occurred while fetching the taxon_name corresponding to id {taxon_id}.')
        return taxon_name
