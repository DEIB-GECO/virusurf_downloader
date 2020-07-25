from Bio import Entrez
from loguru import logger

cached_taxon_id = {}


def host_taxon_id_from_ncbi_taxon_name(taxon_name: str):
    if not taxon_name:
        return None
    else:
        global cached_taxon_id
        taxon_id = cached_taxon_id.get(taxon_name)
        if taxon_id == -1:  # -1 means the cached taxon_id for this taxon name was searched before
            return None
        elif taxon_id is None:
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
        return taxon_id
