import re
from datetime import datetime

from typing import Tuple, Callable, Generator, Optional

from Bio import Entrez
from dateutil.parser import parse
from loguru import logger

import cleaning_module
from data_sources.common_methods_host_sample import host_taxon_id_from_ncbi_taxon_name
from data_sources.virus_sample import VirusSample

gene_protein_name_replacements = {
    "E": ("E", "E (envelope protein)"),
    "M": ("M", "M (membrane glycoprotein)"),
    "N": ("N", "N (nucleocapsid phosphoprotein)"),
    "NSP16": ("ORF1ab", "NSP16 (2'-O-ribose methyltransferase)"),
    "NSP5": ("ORF1ab", "NSP5 (3C-like proteinase)"),
    "NSP14": ("ORF1ab", "NSP14 (3'-to-5' exonuclease)"),
    "NSP15": ("ORF1ab", "NSP15 (endoRNAse)"),
    "NSP13": ("ORF1ab", "NSP13 (helicase)"),
    "NSP1": ("ORF1ab", "NSP1 (leader protein)"),
    "NSP10": ("ORF1ab", "NSP10"),
    "NSP11": ("ORF1ab", "NSP11"),
    "NSP2": ("ORF1ab", "NSP2"),
    "NSP3": ("ORF1ab", "NSP3"),
    "NSP4": ("ORF1ab", "NSP4"),
    "NSP6": ("ORF1ab", "NSP6"),
    "NSP7": ("ORF1ab", "NSP7"),
    "NSP8": ("ORF1ab", "NSP8"),
    "NSP9": ("ORF1ab", "NSP9"),
    "NSP12": ("ORF1ab", "NSP12 (RNA-dependent RNA polymerase)"),
    "NS3": ("ORF3a", "NS3 (ORF3a protein)"),
    "NS6": ("ORF6", "NS6 (ORF6 protein)"),
    "NS7a": ("ORF7a", "NS7a (ORF7a protein)"),
    "NS7b": ("ORF7b", "NS7b (ORF7b)"),
    "NS8": ("ORF8", "NS8 (ORF8 protein)"),
    "Spike": ("S", "Spike (surface glycoprotein)")
}


def replace_gene_and_protein_name(protein_name):
    return gene_protein_name_replacements.get(protein_name, (None, protein_name))


class GISAIDSarsCov2Sample(VirusSample):
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)
    aa_variant_regex = re.compile(r'(\D+)(\d+)(\D+)')

    def __init__(self, sequence_dict):
        super().__init__()
        self.sequence_dict = sequence_dict

    def internal_id(self):
        return self.primary_accession_number() or str(self.sequence_dict)

    def primary_accession_number(self):
        try:
            return self.sequence_dict['covv_accession_id']
        except KeyError:
            return None

    def alternative_accession_number(self):
        return None

    def strain(self):
        try:
            return self.sequence_dict['covv_virus_name']
        except KeyError:
            return None

    def is_reference(self):
        try:
            return self.sequence_dict['is_reference']
        except KeyError:
            return None

    def is_complete(self):
        try:
            return self.sequence_dict['is_complete']
        except KeyError:
            return None

    def nucleotide_sequence(self):
        return None

    def strand(self):
        try:
            return strip_or_none(self.sequence_dict['covv_strand'])
        except KeyError:
            return None

    def length(self):
        try:
            return self.sequence_dict['sequence_length']
        except KeyError:
            return None

    def gc_percent(self):
        try:
            return self.sequence_dict['gc_content']
        except KeyError:
            return None

    def n_percent(self):
        return self.sequence_dict.get('n_content')

    def lineage(self):
        try:
            return strip_or_none(self.sequence_dict['covv_lineage'])
        except KeyError:
            return None

    def clade(self):
        try:
            clade = self.sequence_dict['covv_clade']
            return strip_or_none(clade)
        except KeyError:
            return None

    def sequencing_technology(self):
        return None

    def assembly_method(self):
        return None

    def coverage(self):
        return None

    def collection_date(self):
        collection_date = self.sequence_dict.get('covv_collection_date')
        if collection_date:
            collection_date = parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
        return collection_date

    def isolation_source(self):
        try:
            isol_source = self.sequence_dict['covv_specimen'].strip().lower()
            if isol_source in ['unknown', 'other', '']:
                return None
            else:
                return strip_or_none(isol_source)
        except KeyError:
            return None
        except AttributeError:
            return None

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        country, region, geo_group = None, None, None
        try:
            covv_location = self.sequence_dict['covv_location']
            split_locations = [x.strip() for x in covv_location.split('/')]
            geo_group = split_locations[0]
            country = split_locations[1]
            region = split_locations[2]
        except KeyError:
            pass
        except IndexError:
            pass
        return strip_or_none(country), strip_or_none(region), strip_or_none(geo_group)

    def submission_date(self):
        submission_date = self.sequence_dict.get('covv_subm_date')
        if submission_date:
            submission_date = parse(submission_date, default=self.default_datetime).strftime('%Y-%m-%d')
        return submission_date

    def originating_lab(self) -> Optional[str]:
        lab = self.sequence_dict.get('covv_orig_lab')
        return strip_or_none(lab)

    def sequencing_lab(self) -> Optional[str]:
        lab = self.sequence_dict.get('covv_subm_lab')
        return strip_or_none(lab)

    def host_taxon_name(self) -> Optional[str]:
        taxon_name = self.sequence_dict.get('covv_host')
        if taxon_name is not None:
            taxon_name = taxon_name.strip()
            taxon_name = cleaning_module.correct_typos(taxon_name)
        return taxon_name

    def host_taxon_id(self) -> Optional[int]:
        return host_taxon_id_from_ncbi_taxon_name(self.host_taxon_name())

    def gender(self) -> Optional[str]:
        return None

    def age(self) -> Optional[int]:
        return None

    def database_source(self) -> str:
        return 'GISAID'

    def bioproject_id(self):
        return None

    def nucleotide_variants_and_effects(self, aligner: Callable) -> Generator[Tuple, None, None]:
        yield from ()   # empty generator

    def annotations_and_amino_acid_variants(self, reference_virus_sample=None) -> Generator[Tuple, None, None]:
        start, stop, feature_type, gene_name, db_xref_merged, amino_acid_sequence = None, None, 'CDS', None, None, None
        try:
            mutations = self.sequence_dict['covsurver_prot_mutations']
            if mutations:
                mutations = mutations \
                    .replace('(', '') \
                    .replace(')', '') \
                    .split(',')
                mutations = [string for string in mutations if string and string.lower() != "null"]  # drop empty strings
                # these two are going to contain the result
                annotations, parsed_mutations = [], []
                for mut in mutations:
                    split_fields = mut.split('_')
                    assert len(split_fields) == 2, f'The amino acid mutations can contain more than one _ each! : {mut}'
                    product, variant = split_fields
                    variant_search = GISAIDSarsCov2Sample.aa_variant_regex.search(variant)
                    if variant_search:
                        original_aa = variant_search.group(1)
                        start_pos = variant_search.group(2)
                        alternative_aa = variant_search.group(3)

                        if original_aa.lower() == 'ins':
                            original_aa = '-'
                            _type = 'INS'
                        elif alternative_aa.lower() == 'del':
                            alternative_aa = '-'
                            _type = 'DEL'
                        else:
                            _type = 'SUB'

                        length = max(len(original_aa), len(alternative_aa))

                        # if annotation with this product doesn't exists, add it
                        try:
                            product_idx = annotations.index(product)
                        except ValueError:
                            annotations.append(product)
                            product_idx = len(annotations) - 1
                        # add mutation to the corresponding group of annotation products
                        try:
                            already_present_mutations = parsed_mutations[product_idx]
                        except IndexError:
                            already_present_mutations = []
                            parsed_mutations.insert(product_idx, already_present_mutations)
                        already_present_mutations.append((original_aa, alternative_aa, start_pos, length, _type))

                for product, _mutations in zip(annotations, parsed_mutations):
                    gene_name, product = replace_gene_and_protein_name(product)
                    yield start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, _mutations
        except KeyError:
            yield from ()   # empty generator

    def on_before_multiprocessing(self):
        pass

    def nucleotide_var_aligner(self):
        return lambda: None


def strip_or_none(string_or_none: Optional[str]):
    if string_or_none is not None:
        return string_or_none.strip()
    else:
        return None