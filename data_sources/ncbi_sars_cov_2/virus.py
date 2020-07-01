import os
from collections import OrderedDict

from Bio import Entrez
from loguru import logger
from lxml import etree
from lxml.etree import XMLSyntaxError
from tqdm import tqdm
from typing import List

from data_sources.ncbi_sars_cov_2.sample import NCBISarsCov2Sample
from locations import get_local_folder_for, FileType

from xml_helper import text_at_node
from data_sources.virus import VirusSource


# noinspection PyMethodMayBeStatic
class NCBISarsCov2(VirusSource):

    name = 'NCBI_sars_cov_2'

    def __init__(self):
        super().__init__()
        logger.info(f'importing virus {self.name}')
        self.tax_tree = download_virus_taxonomy_as_xml(
            get_local_folder_for(source_name=self.name, _type=FileType.TaxonomyData),
            self.taxon_id())

    def taxon_id(self):
        # taxon_id = int(text_at_node(self.tax_tree, './Taxon/TaxId'))
        return 2697049

    def taxon_name(self):
        return text_at_node(self.tax_tree, './Taxon/ScientificName')

    def family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "family"]/ScientificName')

    def sub_family(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "subfamily"]/ScientificName')

    def genus(self):
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "genus"]/ScientificName')

    def species(self):
        # species_taxon_id = text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/TaxId')
        return text_at_node(self.tax_tree, './/LineageEx/Taxon[./Rank/text() = "species"]/ScientificName', mandatory=False)

    def equivalent_names(self):
        genbank_acronym = text_at_node(self.tax_tree, './/GenbankAcronym', mandatory=False)
        equivalent_names = self.tax_tree.xpath('.//EquivalentName')
        equivalent_names = [x.text for x in equivalent_names]
        if genbank_acronym:
            equivalent_names.insert(0, genbank_acronym)
        equivalent_names = list(OrderedDict.fromkeys(equivalent_names))
        equivalent_names = ", ".join(equivalent_names)
        return equivalent_names

    def molecule_type(self):
        return 'RNA'

    def is_single_stranded(self):
        return True

    def is_positive_stranded(self):
        return True

    def virus_samples(self):
        logger.info(f'getting accession ids for virus sequences')
        refseq_accession_id, non_refseq_accession_ids = get_virus_sample_accession_ids(self.taxon_id())
        # logger.warning('Sequence accession ids are hardcoded in main.py.')
        # refseq_accession_id = 1798174254      # hardcoded value for tests
        # non_refseq_accession_ids = [1852393386, 1852393360, 1852393373,  1852393398, 1859035944, 1859035892, 1800242649, 1800242657, 1858732896, 1799706760, 1800242651, 1800242659, 1858732909, 1799706762, 1800242653, 1800242661, 1858732922, 1800242639, 1800242655, 1850952215, 1859094271]
        # non_refseq_accession_ids = non_refseq_accession_ids[0:2]      # limit number of seq
        # non_refseq_accession_ids = non_refseq_accession_ids[::-1]     # invert seq
        sequence_accession_ids = [refseq_accession_id] + non_refseq_accession_ids

        sample_local_download_dir = get_local_folder_for(source_name=self.name, _type=FileType.SequenceOrSampleData)

        for virus_seq_acc_id in tqdm(sequence_accession_ids):
            try:
                virus_sample_file_path = download_or_get_virus_sample_as_xml(sample_local_download_dir, virus_seq_acc_id)
                yield NCBISarsCov2Sample(virus_sample_file_path, virus_seq_acc_id)
            except XMLSyntaxError:
                logger.error(f'virus sample file with accession id {virus_seq_acc_id} was malformed or empty and it was deleted.')
                delete_virus_sample_xml(sample_local_download_dir, virus_seq_acc_id)


def download_virus_taxonomy_as_xml(containing_directory: str, taxon_id) -> etree.ElementTree:
    # write taxonomy tree
    local_file_path = f"{containing_directory}{os.path.sep}{taxon_id}.xml"
    if not os.path.exists(local_file_path):
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype=None, retmode="xml") as handle, \
                open(local_file_path, 'w') as f:
            for line in handle:
                f.write(line)
    tax_tree = etree.parse(local_file_path, parser=etree.XMLParser(remove_blank_text=True))
    return tax_tree


def get_virus_sample_accession_ids(virus_specie_taxon_id: int) -> (int, List[int]):
    """
    For a virus specie, it returns the accession ids of the reference sample (refseq) and a list of ids for the non-reference
    samples.
    """
    # get reference sample accession id
    handle = Entrez.esearch(db="nuccore", term=f"(txid{virus_specie_taxon_id}[Organism]) AND srcdb_refseq[Properties]")
    response = Entrez.read(handle)
    handle.close()
    # Example response
    # {
    #     "Count": "1",                       <-- number of total records matching the query
    #     "RetMax": "1",
    #     "RetStart": "0",
    #     "IdList": ["1798174254"],            <-- accession id of refseq
    #     "TranslationSet": [],
    #     "TranslationStack": [
    #         {"Term": "txid2697049[Organism]", "Field": "Organism", "Count": "5511", "Explode": "Y"},
    #         {"Term": "srcdb_refseq[Properties]", "Field": "Properties", "Count": "66995641", "Explode": "N"},
    #         "AND"
    #     ],
    #     "QueryTranslation": "txid2697049[Organism] AND srcdb_refseq[Properties]"
    # }
    assert int(response['Count']) == 1, \
        "no reference sample found or multiple RefSeqs" + "please check: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch"
    refseq_accession_id = response['IdList'][0]

    # get non-refseq accession ids
    #   total number of sequences
    with Entrez.esearch(db="nuccore",
                        term=f"(txid{virus_specie_taxon_id}[Organism]) NOT srcdb_refseq[Properties]",
                        rettype='count') as handle:
        response = Entrez.read(handle)
        total_records = int(response['Count'])
    #   do pagination
    non_refseq_accessions_ids = list()
    RECORDS_PER_PAGE = 1000
    page_number = 0
    with tqdm(total=total_records) as progress_bar:
        while total_records > page_number*RECORDS_PER_PAGE:
            with Entrez.esearch(db="nuccore", term=f"(txid{virus_specie_taxon_id}[Organism]) NOT srcdb_refseq[Properties]",
                                retmax=RECORDS_PER_PAGE, retstart=page_number*RECORDS_PER_PAGE) as handle:
                response = Entrez.read(handle)
            for x in response['IdList']:
                non_refseq_accessions_ids.append(int(x))
            page_number += 1
            progress_bar.update(RECORDS_PER_PAGE if page_number*RECORDS_PER_PAGE < total_records else total_records-((page_number-1)*RECORDS_PER_PAGE))
    assert len(non_refseq_accessions_ids) == total_records, 'Some of the non-refseq accession ids were not correctly downloaded'
    return refseq_accession_id, non_refseq_accessions_ids


def download_or_get_virus_sample_as_xml(containing_directory: str, sample_accession_id: int) -> str:
    """
    :param containing_directory: directory where the file will be downloaded and cached
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
    with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml") as handle, open(local_file_path, 'w') as f:
        for line in handle:
            f.write(line)
    return local_file_path


def delete_virus_sample_xml(containing_directory: str, sample_accession_id: int):
    """
    :param containing_directory: directory where the file resides
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    """
    local_file_path = f"{containing_directory}{os.path.sep}{sample_accession_id}.xml"
    try:
        os.remove(local_file_path)
    except OSError as e:
        logger.error(f"Failed to remove file {local_file_path} with error: {e.strerror}")
