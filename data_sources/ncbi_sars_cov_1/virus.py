import os
from lxml.etree import XMLSyntaxError

from typing import List

from Bio import Entrez
from loguru import logger
from tqdm import tqdm

from data_sources.ncbi_sars_cov_1.sample import NCBISarsCov1Sample
from data_sources.ncbi_sars_cov_2.virus import NCBISarsCov2, download_or_get_virus_sample_as_xml, \
    delete_virus_sample_xml
from locations import get_local_folder_for, FileType


class NCBISarsCov1(NCBISarsCov2):

    name = 'NCBI_sars_cov_1'

    def __init__(self):
        super().__init__()

    def taxon_id(self):
        return 694009

    def taxon_name(self):
        return super().taxon_name()

    def family(self):
        return super().family()

    def sub_family(self):
        return super().sub_family()

    def genus(self):
        return super().genus()

    def species(self):
        return self.taxon_name()

    def equivalent_names(self):
        return super().equivalent_names()

    def molecule_type(self):
        return super().molecule_type()

    def is_single_stranded(self):
        return super().is_single_stranded()

    def is_positive_stranded(self):
        return super().is_positive_stranded()

    def virus_samples(self):
        refseq_accession_id, non_refseq_accession_ids = get_virus_sample_accession_ids(self.taxon_id())
        # refseq_accession_id = 30271926
        # non_refseq_accession_ids = [1773397070, 1773397068, 1773397066, 1773397050, 1773397048, 1773397046, 1773397044, 1773397042, 1773397024, 1773397022, 1773397020, 1773397018, 1773397002, 1773397000, 1773396998, 1773396996, 1773396994, 1773396992]
        non_refseq_accession_ids = non_refseq_accession_ids[:20]
        sequence_accession_ids = [refseq_accession_id] + non_refseq_accession_ids

        sample_local_download_dir = get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)

        for virus_seq_acc_id in tqdm(sequence_accession_ids):
            try:
                virus_sample_file_path = download_or_get_virus_sample_as_xml(sample_local_download_dir, virus_seq_acc_id)
                yield NCBISarsCov1Sample(virus_sample_file_path, virus_seq_acc_id)
            except XMLSyntaxError:
                logger.error(f'virus sample file with accession id {virus_seq_acc_id} was malformed or empty and it was deleted.')
                delete_virus_sample_xml(sample_local_download_dir, virus_seq_acc_id)


def get_virus_sample_accession_ids(virus_specie_taxon_id: int) -> (int, List[int]):
    """
    For a virus specie, it returns the accession ids of the reference sample (refseq) and a list of ids for the non-reference
    samples.
    """
    # get reference sample accession id
    handle = Entrez.esearch(db="nuccore", term=f'txid694009[Organism] AND srcdb_refseq[Properties] NOT txid2697049[Organism]')
    # the query is similar to the one for SarsCov2 + the addition of the filter Publication Date <= 2018 to exclude the Sars-Cov2 reference sequence from the result list.
    response = Entrez.read(handle)
    handle.close()
    # Example response
    # {
    #     "Count": "1",                             <-- number of total records matching the query
    #     "RetMax": "1",
    #     "RetStart": "0",
    #     "IdList": ["30271926"],                   <-- accession id of refseq
    #     "TranslationSet": [],
    #     "TranslationStack": [
    #         {"Term": "txid694009[Organism]", "Field": "Organism", "Count": "9640", "Explode": "Y"},
    #         {"Term": "srcdb_refseq[Properties]", "Field": "Properties", "Count": "67290614", "Explode": "N" },
    #         "AND",
    #         {"Term": "0001/01/01\"[PDAT]", "Field": "PDAT", "Count": "0", "Explode": "N" },
    #         {"Term": "2018/12/31\"[PDAT]", "Field": "PDAT", "Count": "0", "Explode": "N" },
    #         "RANGE",
    #         "AND"
    #     ],
    #     "QueryTranslation": "txid694009[Organism] AND srcdb_refseq[Properties] AND \"0001/01/01\"[PDAT] : \"2018/12/31\"[PDAT]"
    # }
    assert int(response['Count']) == 1, \
        "no reference sample found or multiple RefSeqs" + "please check: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch"
    refseq_accession_id = response['IdList'][0]

    # get non-refseq accession ids
    #   total number of sequences
    with Entrez.esearch(db="nuccore",
                        term=f'txid694009[Organism:noexp] NOT srcdb_refseq[Properties] NOT txid2697049[Organism]',
                        rettype='count') as handle:
        response = Entrez.read(handle)
        total_records = int(response['Count'])
    #   do pagination
    non_refseq_accessions_ids = list()
    RECORDS_PER_PAGE = 1000
    page_number = 0
    with tqdm(total=total_records) as progress_bar:
        while total_records > page_number*RECORDS_PER_PAGE:
            with Entrez.esearch(db="nuccore", term=f"txid694009[Organism:noexp] NOT srcdb_refseq[Properties] NOT txid2697049[Organism]",
                                retmax=RECORDS_PER_PAGE, retstart=page_number*RECORDS_PER_PAGE) as handle:
                response = Entrez.read(handle)
            for x in response['IdList']:
                non_refseq_accessions_ids.append(int(x))
            page_number += 1
            progress_bar.update(RECORDS_PER_PAGE if page_number*RECORDS_PER_PAGE < total_records else total_records-((page_number-1)*RECORDS_PER_PAGE))
    assert len(non_refseq_accessions_ids) == total_records, 'Some of the non-refseq accession ids were not correctly downloaded'
    return refseq_accession_id, non_refseq_accessions_ids
