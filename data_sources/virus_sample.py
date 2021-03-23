from typing import Tuple, Callable, Generator, Optional


class VirusSample:

    def __init__(self):
        pass

    def internal_id(self):
        """
        An id used internally by this program when needed to log some errors. Choose whatever form and type you like.
        This id is not going to be written in the output data and possibly it will be printed only on the console in order
        to trace errors.
        """
        raise NotImplementedError

    def primary_accession_number(self):
        raise NotImplementedError

    def alternative_accession_number(self):
        raise NotImplementedError

    def strain(self):
        raise NotImplementedError

    def is_reference(self):
        raise NotImplementedError

    def is_complete(self):
        raise NotImplementedError

    def nucleotide_sequence(self):
        raise NotImplementedError

    def strand(self):
        raise NotImplementedError

    def length(self):
        raise NotImplementedError

    def gc_percent(self):
        raise NotImplementedError

    def n_percent(self):
        raise NotImplementedError

    def lineage(self):
        raise NotImplementedError

    def clade(self):
        raise NotImplementedError

    def sequencing_technology(self):
        raise NotImplementedError

    def assembly_method(self):
        raise NotImplementedError

    def coverage(self):
        raise NotImplementedError

    def collection_date(self):
        raise NotImplementedError

    def isolation_source(self):
        raise NotImplementedError

    def province__region__country__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        raise NotImplementedError

    def submission_date(self):
        raise NotImplementedError

    def originating_lab(self) -> Optional[str]:
        raise NotImplementedError

    def sequencing_lab(self) -> Optional[str]:
        raise NotImplementedError

    def host_taxon_name(self) -> Optional[str]:
        raise NotImplementedError

    def host_taxon_id(self) -> Optional[int]:
        raise NotImplementedError

    def gender(self) -> Optional[str]:
        raise NotImplementedError

    def age(self) -> Optional[int]:
        raise NotImplementedError

    def database_source(self) -> str:
        raise NotImplementedError

    def bioproject_id(self):
        raise NotImplementedError

    def nucleotide_variants_and_effects(self, aligner: Callable) -> Generator[Tuple, None, None]:
        raise NotImplementedError

    def annotations_and_amino_acid_variants(self, reference_virus_sample) -> Generator[Tuple, None, None]:
        raise NotImplementedError

    def on_before_multiprocessing(self):
        raise NotImplementedError

    def nucleotide_var_aligner(self):
        raise NotImplementedError


# class DoNotImportSample(RollbackTransactionAndRaise):
#     def __init__(self, sample_name_or_accession_id):
#         self.internal_name_or_accession_id = sample_name_or_accession_id
#         super().__init__(f'Request to abort the import of sample with accession id {sample_name_or_accession_id}')
