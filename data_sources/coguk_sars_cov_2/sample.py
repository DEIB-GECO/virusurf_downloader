import os
from loguru import logger
from typing import Generator, Tuple, Callable, Optional, List, Iterable

from data_sources.virus_sample import VirusSample
from locations import local_folder_nuc_variant_and_effects


class COGUKSarsCov2Sample(VirusSample):
    # DICT KEYS
    STRAIN_NAME = 'strain_name'
    NUC_SEQUENCE = 'nuc_sequence'
    METADATA_RAW_STRING = 'meta_string'
    COUNTRY = 'country'
    REGION = 'region'
    COLLECTION_DATE = 'collection_date'
    LINEAGE = 'lineage'

    def __init__(self, sample_dict):
        super().__init__()
        self.sample_dict = self.parse_sample_dict(sample_dict)

    def parse_sample_dict(self, sample_dict):
        meta: str = sample_dict.get(self.METADATA_RAW_STRING)
        if meta:
            try:
                country, region, collection_date, epi_week, lineage, lineage_support = meta.split(',')
            except:
                logger.exception(f'Unable to parse metadata for sample {sample_dict[self.STRAIN_NAME]}. Metadata string is:\n\t{meta}')
                return sample_dict
            else:
                output = {
                    self.STRAIN_NAME: sample_dict['strain_name'],
                    self.NUC_SEQUENCE: sample_dict['nuc_sequence'],
                    self.COUNTRY: country,
                    self.REGION: region,
                    self.COLLECTION_DATE: collection_date,
                    self.LINEAGE: lineage
                }
                # forget about epi_week and lineage_support
                return output
        else:
            return sample_dict

    def internal_id(self):
        return self.strain().replace('/', '-')

    def primary_accession_number(self):
        return self.strain()

    def alternative_accession_number(self):
        return None

    def strain(self):
        return self.sample_dict[self.STRAIN_NAME]

    def is_reference(self):
        return False

    def is_complete(self):
        return None

    def nucleotide_sequence(self):
        return self.sample_dict[self.NUC_SEQUENCE]

    def strand(self):
        return None

    def length(self):
        return None

    def gc_percent(self):
        return None

    def n_percent(self):
        return None

    def lineage(self):
        return self.sample_dict.get(self.LINEAGE)

    def clade(self):
        return None

    def sequencing_technology(self):
        return None

    def assembly_method(self):
        return None

    def coverage(self):
        return None

    def collection_date(self):
        return self.sample_dict.get(self.COLLECTION_DATE)

    def isolation_source(self):
        return None

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        return self.sample_dict.get(self.COUNTRY), self.sample_dict.get(self.REGION), None

    def submission_date(self):
        return None

    def originating_lab(self) -> Optional[str]:
        return None

    def sequencing_lab(self) -> Optional[str]:
        return None

    def taxon_name(self) -> Optional[str]:
        return 'Human'

    def taxon_id(self) -> Optional[int]:
        return 9606

    def gender(self) -> Optional[str]:
        return None

    def age(self) -> Optional[int]:
        return None

    def database_source(self) -> str:
        return 'COG-UK'

    def bioproject_id(self):
        return None

    def _call_or_read_nuc_variants_and_effects(self, aligner: Callable) -> Generator[Iterable, None, None]:
        cache_file_path = f'{local_folder_nuc_variant_and_effects}{os.path.sep}{self.internal_id()}'
        try:
            with open(cache_file_path, mode='r') as f:
                logger.trace(f'reading nucleotide variants for sample {self.internal_id()} from disk')
                for line in f:
                    yield line.rstrip().split('\t')
        except FileNotFoundError:
            # compute them and write to file
            logger.trace(f'calling nucleotide variants for sample {self.internal_id()}...')
            variants: Tuple = aligner(sequence=self.nucleotide_sequence(), sequence_id=self.internal_id())
            genbank_annotated_variants = variants[0]
            with open(cache_file_path, mode='w') as f:
                for variant in genbank_annotated_variants:
                    _, start_original, _, _, _, _, others, snpeff_ann = variant.split("\t")
                    f.write(f'{start_original}\t{others}\t{snpeff_ann}')  # snpeff_ann already contains a \n
                    yield [start_original, others, snpeff_ann]

    def nucleotide_variants_and_effects(self, aligner: Callable) -> Generator[Tuple, None, None]:
        """
        Return a generator that iterates over the variants produced by snpEff and each time it yields a tuple of values
        describing a row of the table NucleotideVariant. The tuple reports in order:
        sequence original, sequence alternative, start original, start alternative, variant length, variant type.
        """
        for start_original, others, snpeff_ann in self._call_or_read_nuc_variants_and_effects(aligner):
            # get a variant
            variant_type, start_alternative, variant_length, sequence_original, sequence_alternative = others.split(',')
            # get variant effects and remove duplicated effects
            variant_impacts = set()
            for ann in snpeff_ann.split(","):
                s = ann.split("|")
                variant_impacts.add((s[1], s[2], s[3]))
            # return each variant with associated effects
            yield sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type, variant_impacts

    def annotations_and_amino_acid_variants(self, reference_virus_sample) -> Generator[Tuple, None, None]:
        yield from ()   # empty generator

    def on_before_multiprocessing(self):
        pass

    def nucleotide_var_aligner(self):
        return lambda: None
