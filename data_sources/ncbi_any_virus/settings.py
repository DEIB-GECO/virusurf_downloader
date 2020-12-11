from typing import Optional
from os.path import sep


reference_sample_query: Optional[str] = None # initialized elsewhere
non_reference_samples_query: Optional[str] = None # initialized elsewhere
virus_taxon_id: Optional[int] = None # initialized elsewhere
chromosome_name: Optional[str] = None # initialized elsewhere
snpeff_db_name: Optional[str] = None # initialized elsewhere
annotation_file_path: Optional[str] = None # initialized elsewhere
log_with_name: Optional[str] = None # initialized elsewhere
generated_dir_name: Optional[str] = None # initialized elsewhere

nucleotide_reference_sequence: Optional[str] = None  # initialized elsewhere
length_of_nuc_reference_sequence: Optional[int] = None  # initialized elsewhere

def initialize_settings(_reference_sample_query: str, _non_reference_samples_query: str, _virus_taxon_id: int,
                        _chromosome_name: str, _snpeff_db_name: str, _annotation_file_path: str, _log_with_name: str,
                        _generated_dir_name: str):
    global reference_sample_query, non_reference_samples_query, virus_taxon_id, chromosome_name, snpeff_db_name, \
        annotation_file_path, log_with_name, generated_dir_name

    reference_sample_query = _reference_sample_query
    non_reference_samples_query = _non_reference_samples_query
    virus_taxon_id = _virus_taxon_id
    chromosome_name = _chromosome_name
    snpeff_db_name = _snpeff_db_name
    annotation_file_path = _annotation_file_path
    log_with_name = _log_with_name
    generated_dir_name = _generated_dir_name


def set_nucleotide_refseq(seq: str):
    global nucleotide_reference_sequence, length_of_nuc_reference_sequence
    nucleotide_reference_sequence = seq
    length_of_nuc_reference_sequence = len(seq)


known_settings = {
    "sars_cov_2": {
        "reference_sample_query": "txid2697049[Organism] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid2697049[Organism] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 2697049,
        "chromosome_name": "NC_045512",
        "snpeff_db_name": "new_ncbi_sars_cov_2",
        "annotation_file_path": f".{sep}annotations{sep}new_ncbi_sars_cov_2.tsv",
        "log_with_name": "NCBI-SC2",
        "generated_dir_name": "ncbi_sars_cov_2"
    },
    "sars_cov_1": {
        "reference_sample_query": "txid694009[Organism] NOT txid2697049[Organism] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid694009[Organism:noexp] NOT txid2697049[Organism]",
        "virus_taxon_id": 694009,
        "chromosome_name": "NC_004718",
        "snpeff_db_name": "sars_cov_1",
        "annotation_file_path": f".{sep}annotations{sep}sars_cov_1.tsv",
        "log_with_name": "SC1",
        "generated_dir_name": "sars_cov_1"
    },
    "dengue_1": {
        "reference_sample_query": "txid11053[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid11053[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 11053,
        "chromosome_name": "NC_001477",
        "snpeff_db_name": "dengue_virus_1",
        "annotation_file_path": f".{sep}annotations{sep}dengue_virus_1.tsv",
        "log_with_name": "DENGUE-1",
        "generated_dir_name": "dengue_1"
    },
    "dengue_2": {
        "reference_sample_query": "txid11060[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid11060[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 11060,
        "chromosome_name": "NC_001474",
        "snpeff_db_name": "dengue_virus_2",
        "annotation_file_path": f".{sep}annotations{sep}dengue_virus_2.tsv",
        "log_with_name": "DENGUE-2",
        "generated_dir_name": "dengue_2"
    },
    "dengue_3": {
        "reference_sample_query": "txid11069[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid11069[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 11069,
        "chromosome_name": "NC_001475",
        "snpeff_db_name": "dengue_virus_3",
        "annotation_file_path": f".{sep}annotations{sep}dengue_virus_3.tsv",
        "log_with_name": "DENGUE-3",
        "generated_dir_name": "dengue_3"
    },
    "dengue_4": {
        "reference_sample_query": "txid11070[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid11070[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 11070,
        "chromosome_name": "NC_002640",
        "snpeff_db_name": "dengue_virus_4",
        "annotation_file_path": f".{sep}annotations{sep}dengue_virus_4.tsv",
        "log_with_name": "DENGUE-4",
        "generated_dir_name": "dengue_4"
    },
    "mers": {
        "reference_sample_query": "txid1335626[Organism:noexp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid1335626[Organism:noexp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 1335626,
        "chromosome_name": "NC_019843",
        "snpeff_db_name": "mers",
        "annotation_file_path": f".{sep}annotations{sep}mers.tsv",
        "log_with_name": "MERS",
        "generated_dir_name": "mers"
    },
    "betacoronavirus_england_1": {
        "reference_sample_query": "txid1263720[Organism:noexp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid1263720[Organism:noexp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 1263720,
        "chromosome_name": "NC_038294",
        "snpeff_db_name": "betacoronavirus_england_1",
        "annotation_file_path": f".{sep}annotations{sep}betacoronavirus_england_1.tsv",
        "log_with_name": "BETA-COV-EN-1",
        "generated_dir_name": "betacov_england_1"
    },
    "zaire_ebolavirus": {
        "reference_sample_query": "txid186538[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid186538[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 186538,
        "chromosome_name": "NC_002549",
        "snpeff_db_name": "zaire_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}zaire_ebolavirus.tsv",
        "log_with_name": "ZAIRE-EV",
        "generated_dir_name": "zaire_ebolavirus"
    },
    "sudan_ebolavirus": {
        "reference_sample_query": "txid186540[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid186540[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 186540,
        "chromosome_name": "NC_006432",
        "snpeff_db_name": "sudan_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}sudan_ebolavirus.tsv",
        "log_with_name": "SUDAN-EV",
        "generated_dir_name": "sudan_ebolavirus"
    },
    "reston_ebolavirus": {
        "reference_sample_query": "txid186539[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid186539[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 186539,
        "chromosome_name": "NC_004161",
        "snpeff_db_name": "reston_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}reston_ebolavirus.tsv",
        "log_with_name": "RESTON-EV",
        "generated_dir_name": "reston_ebolavirus"
    },
    "bundibugyo_ebolavirus": {
        "reference_sample_query": "txid565995[Organism:noexp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid565995[Organism:noexp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 565995,
        "chromosome_name": "NC_014373",
        "snpeff_db_name": "bundibugyo_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}bundibugyo_ebolavirus.tsv",
        "log_with_name": "BUNDIBUGYO_EV",
        "generated_dir_name": "bundibugyo_ebolavirus"
    },
    "bombali_ebolavirus": {
        "reference_sample_query": "txid2010960[Organism:noexp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid2010960[Organism:noexp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 2010960,
        "chromosome_name": "NC_039345",
        "snpeff_db_name": "bombali_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}bombali_ebolavirus.tsv",
        "log_with_name": "BOMBALI-EV",
        "generated_dir_name": "bombali_ebolavirus"
    },
    "tai_forest_ebolavirus": {
        "reference_sample_query": "txid186541[Organism:exp] AND srcdb_refseq[Properties]",
        "non_reference_samples_query": "txid186541[Organism:exp] NOT srcdb_refseq[Properties]",
        "virus_taxon_id": 186541,
        "chromosome_name": "NC_014372",
        "snpeff_db_name": "tai_forest_ebolavirus",
        "annotation_file_path": f".{sep}annotations{sep}tai_forest_ebolavirus.tsv",
        "log_with_name": "TAI-EV",
        "generated_dir_name": "tai_forest_ebolavirus"
    }
}
known_settings["new_ncbi_sars_cov_2"] = known_settings["sars_cov_2"]