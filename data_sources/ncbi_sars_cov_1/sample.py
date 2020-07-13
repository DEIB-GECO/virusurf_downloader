import re
from typing import Generator, Tuple, Callable, Iterable, Optional

from data_sources.ncbi_sars_cov_2.sample import NCBISarsCov2Sample
from xml_helper import text_at_node


class NCBISarsCov1Sample(NCBISarsCov2Sample):
    virus_name = 'NCBI_sars_cov_1'
    cached_taxon_id = {}

    def __init__(self, virus_sample_file_path: str, internal_accession_id):
        super().__init__(virus_sample_file_path, internal_accession_id)

    def internal_id(self):
        return super().internal_id()

    def primary_accession_number(self):
        return super().primary_accession_number()

    def alternative_accession_number(self):
        return super().alternative_accession_number()

    def strain(self):
        return super().strain()

    def is_reference(self):
        return super().is_reference()

    def is_complete(self):
        definition = text_at_node(self.sample_xml, ".//INSDSeq_definition")
        definition_0 = definition.split(";")[0]
        definition_0_last = definition_0.split(",")[-1]
        definition_0_last = definition_0_last.strip()
        if definition_0_last in ['complete genome', ]:
            return True
        elif definition_0_last in ['partial cds', 'complete cds', 'partial genome']:
            return False
        else:
            length = self.length()
            if length and length < 28374:      # 95 % of the length of the reference sequence
                return False
            else:
                return None

    def nucleotide_sequence(self):
        return super().nucleotide_sequence()

    def strand(self):
        return super().strand()

    def length(self):
        return super().length()

    def gc_percent(self):
        return super().gc_percent()

    def n_percent(self):
        return super().n_percent()

    def lineage(self):
        super().lineage()

    def clade(self):
        super().clade()

    def sequencing_technology(self):
        return super().sequencing_technology()

    def assembly_method(self):
        return super().assembly_method()

    def coverage(self):
        return super().coverage()

    def collection_date(self):
        return super().collection_date()

    def isolation_source(self):
        return super().isolation_source()

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        return super().country__region__geo_group()

    def _init_and_get_journal(self):
        if not self._journal:
            try:
                return super()._init_and_get_journal()
            except AssertionError:
                # try to get patent information
                definition = text_at_node(self.sample_xml, './/INSDSeq_definition', mandatory=False)
                if definition and 'patent' in definition.lower():
                    references = self.sample_xml.xpath('.//INSDReference/INSDReference_journal')
                    assert len(references) > 0, 'no journal info for patent'
                    for ref in references:
                        journal_string: str = text_at_node(ref, xpath_string='.', mandatory=True)
                        _journal_parts = re.match(r'^.* (\d+-\w+-\d+)(?: (.+))?$', journal_string)
                        if _journal_parts:
                            self._journal = _journal_parts.groups()
                            break
        return self._journal

    def submission_date(self):
        return super().submission_date()

    def originating_lab(self) -> Optional[str]:
        return super().originating_lab()

    def sequencing_lab(self) -> Optional[str]:
        return super().sequencing_lab()

    def host_taxon_name(self) -> Optional[str]:
        return super().host_taxon_name()

    def host_taxon_id(self) -> Optional[int]:
        return super().host_taxon_id()

    def gender(self) -> Optional[str]:
        return super().gender()

    def age(self) -> Optional[int]:
        return super().age()

    def database_source(self) -> str:
        return super().database_source()

    def bioproject_id(self):
        return super().bioproject_id()

    def nucleotide_variants_and_effects(self, aligner: Callable) -> Generator[Tuple, None, None]:
        return super().nucleotide_variants_and_effects(aligner)

    def annotations_and_amino_acid_variants(self, reference_virus_sample) -> Generator[Tuple, None, None]:
        return super().annotations_and_amino_acid_variants(reference_virus_sample)

    def get_cached_reference_CDS_annotations(self):
        return super().get_cached_reference_CDS_annotations()

    def on_before_multiprocessing(self):
        super().on_before_multiprocessing()

    def nucleotide_var_aligner(self):
        return super().nucleotide_var_aligner()
