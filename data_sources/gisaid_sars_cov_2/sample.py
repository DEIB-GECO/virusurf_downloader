import re
import math
from decimal import Decimal, ROUND_UP
from datetime import datetime
from typing import Tuple, Callable, Generator, Optional, Union
from dateutil.parser import parse
import data_cleaning_module
from data_sources.ncbi_services import host_taxon_id_from_ncbi_taxon_name
from data_sources.virus_sample import VirusSample
from geo_groups import geo_groups

gene_protein_name_replacements = {
    "E": ("E", "E (envelope protein)", "YP_009724392.1"),
    "M": ("M", "M (membrane glycoprotein)", "YP_009724393.1"),
    "N": ("N", "N (nucleocapsid phosphoprotein)", "YP_009724397.2"),
    "NSP16": ("ORF1ab", "NSP16 (2'-O-ribose methyltransferase)", "YP_009725311.1"),
    "NSP5": ("ORF1ab", "NSP5 (3C-like proteinase)", "YP_009742612.1"),
    "NSP14": ("ORF1ab", "NSP14 (3'-to-5' exonuclease)", "YP_009725309.1"),
    "NSP15": ("ORF1ab", "NSP15 (endoRNAse)", "YP_009725310.1"),
    "NSP13": ("ORF1ab", "NSP13 (helicase)", "YP_009725308.1"),
    "NSP1": ("ORF1ab", "NSP1 (leader protein)", "YP_009742608.1"),
    "NSP10": ("ORF1ab", "NSP10", "YP_009742617.1"),
    "NSP11": ("ORF1ab", "NSP11", "YP_009725312.1"),
    "NSP2": ("ORF1ab", "NSP2", "YP_009742609.1"),
    "NSP3": ("ORF1ab", "NSP3", "YP_009742610.1"),
    "NSP4": ("ORF1ab", "NSP4", "YP_009742611.1"),
    "NSP6": ("ORF1ab", "NSP6", "YP_009742613.1"),
    "NSP7": ("ORF1ab", "NSP7", "YP_009742614.1"),
    "NSP8": ("ORF1ab", "NSP8", "YP_009742615.1"),
    "NSP9": ("ORF1ab", "NSP9", "YP_009742616.1"),
    "NSP12": ("ORF1ab", "NSP12 (RNA-dependent RNA polymerase)", "YP_009725307.1"),
    "NS3": ("ORF3a", "NS3 (ORF3a protein)", "YP_009724391.1"),
    "NS6": ("ORF6", "NS6 (ORF6 protein)", "YP_009724394.1"),
    "NS7a": ("ORF7a", "NS7a (ORF7a protein)", "YP_009724395.1"),
    "NS7b": ("ORF7b", "NS7b (ORF7b)", "YP_009725318.1"),
    "NS8": ("ORF8", "NS8 (ORF8 protein)", "YP_009724396.1"),
    "Spike": ("S", "Spike (surface glycoprotein)", "YP_009724390.1")
}


def replace_gene_and_protein_name(protein_name):
    return gene_protein_name_replacements.get(protein_name, (None, protein_name, None))


class GISAIDSarsCov2Sample(VirusSample):
    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)
    aa_variant_regex = re.compile(r'(\D+)(\d+)(\D+)')
    prov_in_parentheses = re.compile(r'.*\((.+)\).*')

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

    def length(self) -> Optional[int]:
        try:
            return int(self.sequence_dict['sequence_length'])
        except KeyError:
            return None

    def gc_percent(self) -> Optional[float]:
        try:
            return round_(self.sequence_dict['gc_content'], 2)
        except KeyError:
            return None

    def n_percent(self) -> Optional[float]:
        try:
            return round_(self.sequence_dict.get('n_content'), 2)
        except KeyError:
            return None

    def lineage(self):
        try:
            lin = strip_or_none(self.sequence_dict['covv_lineage'])
            if lin == 'None':
                return None
            else:
                return lin
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

    def collection_date(self) -> Tuple[Optional[str], Optional[int]]:
        collection_date = self.sequence_dict.get('covv_collection_date')
        if collection_date:
            well_formed_collection_date = parse(collection_date, default=self.default_datetime).strftime('%Y-%m-%d')
            return well_formed_collection_date, collection_date.count('-')
        else:
            return None, None

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

    def province__region__country__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        province, region, country, geo_group = None, None, None, None
        try:
            covv_location = self.sequence_dict['covv_location']
            split_locations = [x.strip() for x in covv_location.split('/')]
            geo_group = split_locations[0]
            country = split_locations[1]
            region = split_locations[2]
            province = split_locations[3]
        except KeyError:
            pass
        except IndexError:
            pass
        #  assign geo group based on the country but use the provided geo_group as fallback value
        geo_group = strip_or_none(geo_group)
        if country:
            country = country.strip()
            geo_group = geo_groups.get(country.lower(), geo_group)
        else:
            country = None

        # clean-up region
        if region:
            region = region.strip()
            if region.endswith('r.'):
                region = region.replace('r.', '')
        else:
            region = None

        # clean-up province
        if province:
            province_has_parentheses = GISAIDSarsCov2Sample.prov_in_parentheses.match(province)
            # if it has parentheses, take what's inside
            if province_has_parentheses:
                province = province_has_parentheses.group(1).strip()
            else:
                # take left part of string if comma is present
                province = province.split(',')[0]
                # replace _ with spaces
                # remove trailing/leading spaces
                province = province.replace('_', ' ').strip().lower()
                province = province.replace('co.', 'county')
            province = province.capitalize()

        return province, region, country, geo_group

    def submission_date(self) -> Optional[str]:
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
            taxon_name = data_cleaning_module.correct_typos(taxon_name)
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
                        start_pos = int(start_pos) if start_pos is not None else None
                        alternative_aa = variant_search.group(3)

                        if original_aa.lower() == 'ins':
                            original_aa = '-'
                            _type = 'INS'
                            alternative_aa = alternative_aa.replace('stop', '*')
                        elif alternative_aa.lower() == 'del':
                            alternative_aa = '-'
                            _type = 'DEL'
                        else:
                            _type = 'SUB'
                            alternative_aa = alternative_aa.replace('stop', '*')

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
                    gene_name, product, db_xref_merged = replace_gene_and_protein_name(product)
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


def round_(number: Union[str, float, int], n_decimals: int) -> float:
    """
    float() are expressed with its closest representation achievable by floating point arithmentic. Sometimes this
    representation introduces strange behaviors, such as when the stored number differs from its string representation,
    or when the rounding is incorrect. To overcome this, python has the Decimal library which instead uses a fixed-point
    arithmetic. Here we use the Decimal class to correctly round a real number, and then we convert it back to float
    in order to restore the expected behavior of a normal equivalence check (this is important when comparing against
    the real values used in a database)
    :param number: the number to round
    :param n_decimals: the decimal precision
    :return: a float correctly rounded
    """
    return float(Decimal(number).quantize(Decimal(str(1 / pow(10, n_decimals))), rounding=ROUND_UP))


def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor
