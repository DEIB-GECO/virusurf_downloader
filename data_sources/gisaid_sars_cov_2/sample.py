import re

from typing import Tuple, Callable, Generator, Optional
from data_sources.virus_sample import VirusSample


class GISAIDSarsCov2Sample(VirusSample):

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
            return self.sequence_dict['covv_strand']
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
        return None

    def lineage(self):
        try:
            return self.sequence_dict['covv_lineage']
        except KeyError:
            return None

    def clade(self):
        try:
            return self.sequence_dict['covv_clade']
        except KeyError:
            return None

    def sequencing_technology(self):
        return None

    def assembly_method(self):
        return None

    def coverage(self):
        return None

    def collection_date(self):
        try:
            return self.sequence_dict['covv_collection_date']
        except KeyError:
            return None

    def isolation_source(self):
        try:
            isol_source = self.sequence_dict['covv_specimen'].strip()
            return isol_source if isol_source else None     # replaces empty string with None
        except KeyError:
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
        return country, region, geo_group

    def submission_date(self):
        try:
            return self.sequence_dict['covv_subm_date']
        except KeyError:
            return None

    def originating_lab(self) -> Optional[str]:
        try:
            return self.sequence_dict['covv_orig_lab']
        except KeyError:
            return None

    def sequencing_lab(self) -> Optional[str]:
        try:
            return self.sequence_dict['covv_subm_lab']
        except KeyError:
            return None

    def taxon_name(self) -> Optional[str]:
        try:
            return self.sequence_dict['covv_host']
        except KeyError:
            return None

    def taxon_id(self) -> Optional[int]:
        return 9606 if 'homo' in self.taxon_name().lower() else None

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

    def annotations_and_amino_acid_variants(self, reference_virus_sample) -> Generator[Tuple, None, None]:
        start, stop, feature_type, gene_name, db_xref_merged, amino_acid_sequence = None, None, 'CDS', None, None, None
        try:
            mutations = self.sequence_dict['covsurver_prot_mutations']
            mutations = mutations \
                .replace('(', '') \
                .replace(')', '') \
                .split(',')
            mutations = [string for string in mutations if string]  # drop empty strings
            # these two are going to contain the result
            annotations, parsed_mutations = [], []
            for mut in mutations:
                split_fields = mut.split('_')
                assert len(split_fields) == 2, f'The amino acid mutations can contain more than one _ each! : {mut}'
                product, variant = split_fields
                variant_search = re.search('(\D+)(\d+)(\D+)', variant)
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
                yield start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence, _mutations
        except KeyError:
            yield from ()   # empty generator

    def on_before_multiprocessing(self):
        pass

    def nucleotide_var_aligner(self):
        return lambda: None
