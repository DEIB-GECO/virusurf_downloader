import re
from collections import Counter
from decimal import Decimal

import lxml
from typing import Tuple, Callable, Generator, List, Iterable
from Bio import Entrez, pairwise2
import numpy as np
from loguru import logger
from lxml import etree
from locations import *
from xml_helper import *
from database_tom import RollbackTransactionWithoutError
from datetime import datetime
from dateutil.parser import parse


def download_or_get_virus_sample_as_xml(sample_accession_id: int) -> str:
    """
    :param sample_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
    :return: the local file path of the download INSDSeq XML file.
    """
    local_file_path = f"{local_folder}/{sample_accession_id}.xml"
    if not os.path.exists(local_file_path):
        with Entrez.efetch(db="nuccore", id=sample_accession_id, rettype="gbc", retmode="xml") as handle, open(local_file_path, 'w') as f:
            for line in handle:
                f.write(line)
    return local_file_path


default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)


class VirusSample:
    """
    Set of getters to take the relevant information from a virus sample in INSDSeq XML format.
    Example of virus sample @ https://www.ncbi.nlm.nih.gov/nuccore/MN908947
    INSDSeq XML format DTD @ https://www.ncbi.nlm.nih.gov/data_specs/dtd/INSD_INSDSeq.mod.dtd
    """
    def __init__(self, virus_sample_file_path: str, internal_accession_id):
        """
        :param virus_sample_file_path: virus sample file (INDSeq XML) path
        :param internal_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
        """
        try:
            self.sample_xml: ElementTree = etree.parse(virus_sample_file_path, parser=etree.XMLParser(remove_blank_text=True))
        except lxml.etree.XMLSyntaxError:
            raise RollbackTransactionWithoutError(f'insertion of sequence with accession id {internal_accession_id} skipped '
                                                  f'because of malformed or empty source file')
        self.internal_accession_id = internal_accession_id

        # init internal variables
        self._journal = None
        self._host_value_list = None
        self._annotations = None


    def underlying_xml_element_tree(self):
        return self.sample_xml

    def primary_accession_number(self):
        return text_at_node(self.sample_xml, './/INSDSeq_accession-version')

    def alternative_accession_number(self):
        return self.internal_accession_id

    def strain(self):
        return text_at_node(self.sample_xml, './/INSDQualifier[./INSDQualifier_name/text() = "strain"]/INSDQualifier_value',
                            False)

    def is_reference(self):
        return text_at_node(self.sample_xml, './/INSDKeyword', mandatory=False) == 'RefSeq'

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
            logger.error(f"In {self.internal_accession_id}, unkown complete string: {definition_0_last}")
            return None

    def nucleotide_sequence(self):
        try:
            return text_at_node(self.sample_xml, './/INSDSeq_sequence')
        except AssertionError:
            raise RollbackTransactionWithoutError(f'Insertion of sequence {self.internal_accession_id} skipped because of missing nucleotide seq.')

    # noinspection PyMethodMayBeStatic
    def strand(self):
        return 'positive'

    def length(self):
        length = int(text_at_node(self.sample_xml, './/INSDSeq_length'))
        assert length == len(self.nucleotide_sequence())
        return length

    def gc_percent(self):
        c = Counter(self.nucleotide_sequence().lower())
        gc_percentage = (c['g'] + c['c']) / (c['g'] + c['c'] + c['a'] + c['t'] + c['u']) * 100
        gc_percentage = Decimal(gc_percentage)
        gc_percentage = round(gc_percentage, 2)
        return gc_percentage

    def n_percent(self):
        c = Counter(self.nucleotide_sequence().lower())
        n_percentage = (c['n']) / self.length() * 100
        n_percentage = Decimal(n_percentage)
        n_percentage = round(n_percentage, 2)
        return n_percentage

    def sequencing_technology(self):
        return structured_comment(self.sample_xml, 'Sequencing Technology')

    def assembly_method(self):
        return structured_comment(self.sample_xml, 'Assembly Method')

    def coverage(self):
        return structured_comment(self.sample_xml, 'Coverage')

    def collection_date(self):
        collection_date = text_at_node(self.sample_xml,
                                       '..//INSDQualifier[./INSDQualifier_name/text() = "collection_date"]/INSDQualifier_value',
                                       mandatory=False)
        if collection_date:
            collection_date = parse(collection_date, default=default_datetime).strftime('%Y-%m-%d')
        return collection_date

    def isolation_source(self):
        source = text_at_node(self.sample_xml,
                              '..//INSDQualifier[./INSDQualifier_name/text() = "isolation_source"]/INSDQualifier_value',
                              mandatory=False)
        if source:
            source = source.lower()
        return source

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        country = None
        region = None
        geo_group = None
        node = text_at_node(self.sample_xml,
                            '..//INSDQualifier[./INSDQualifier_name/text() = "country"]/INSDQualifier_value',
                            mandatory=False)
        if node:
            node = node.split(":")
            country = node[0]
            region = node[1] if len(node) > 1 else None
        return country, region, geo_group

    def _init_and_get_journal(self):
        if not self._journal:
            references = self.sample_xml.xpath('.//INSDReference[./INSDReference_title/text() = "Direct Submission"]')
            assert len(references) > 0, 'there must be at least one direct submission'
            reference = references[0]
            journal_string = text_at_node(reference, "./INSDReference_journal")
            assert journal_string.startswith(
                "Submitted "), 'Cannot find submitted in the Journal of direct submission reference'
            self._journal = re.split("[()]", journal_string, maxsplit=2)
            assert len(self._journal) == 3, f"Journal value problem '{journal_string}' {self._journal}"

            #   # journal details
            #     authors = reference.xpath('.//INSDAuthor')
            #     authors = [x.text for x in authors]
            #     authors = ", ".join(authors)
            #     title = text_at_node(reference, "./INSDReference_title")
            #     journal = text_at_node(reference, "./INSDReference_journal")
            #     publication_date = None
            #     pubmed_id = text_at_node(reference, "./INSDReference_pubmed" , mandatory=False)
            #     popset = None
        return self._journal

    def submitted(self):
        return self._init_and_get_journal()[0]

    def submission_date(self):
        return datetime.strptime(self._init_and_get_journal()[1], '%d-%b-%Y')

    @staticmethod
    def originating_lab() -> Optional[str]:
        return None

    def sequencing_lab(self) -> Optional[str]:
        return self._init_and_get_journal()[2]

    def _init_and_get_host_values(self):
        if not self._host_value_list:
            host_node = text_at_node(self.sample_xml,
                                     './/INSDQualifier[./INSDQualifier_name/text() = "host"]/INSDQualifier_value',
                                     mandatory=False)
            self._host_value_list = [x.strip() for x in host_node.split(";")] if host_node else None
            if self._host_value_list:
                for i in range(1, len(self._host_value_list)):
                    self._host_value_list[i] = self._host_value_list[i].lower()
        return self._host_value_list

    def taxon_name(self) -> Optional[str]:
        host = self._init_and_get_host_values()
        return host[0] if host else None

    def taxon_id(self) -> Optional[int]:
        return 9606 if self.taxon_name() == 'Homo sapiens' else None

    def gender(self) -> Optional[str]:
        host = self._init_and_get_host_values()
        if host:
            gender = 'male' if ('male' in host) else 'female' if 'female' in host else None
            for h in host:
                words = h.split()
                if 'm' in words:
                    gender = 'male'
                    break
                elif 'f' in words:
                    gender = 'female'
                    break
            return gender
        else:
            return None

    def age(self) -> Optional[int]:
        host = self._init_and_get_host_values()
        if host:
            age = next(filter(lambda x: 'age' in x, host), None)
            if age:
                age = int(next(filter(lambda x: x.isdigit(), age.split()), None))
            return age
        else:
            return None

    def database_source(self) -> str:
        return 'RefSeq' if self.is_reference() else 'GenBank'

    def bioproject_id(self):
        return text_at_node(self.sample_xml, './/INSDXref[./INSDXref_dbname/text() = "BioProject"]/INSDXref_id', mandatory=False)

    def _create_or_read_aligned_variants(self, aligner: Callable) -> Generator[Iterable, None, None]:
        cache_file_path = f'{local_folder_variant}{os.path.sep}{self.internal_accession_id}'
        try:
            with open(cache_file_path, mode='r') as f:
                for line in f:
                    yield line.rstrip().split('\t')
        except FileNotFoundError:
            # compute them and write to file
            variants: Tuple = aligner(sequence=self.nucleotide_sequence(), sequence_id=self.internal_accession_id)
            genbank_annotated_variants = variants[0]
            with open(cache_file_path, mode='w') as f:
                for variant in genbank_annotated_variants:
                    _, start_original, _, _, _, _, others, snpeff_ann = variant.split("\t")
                    f.write(f'{start_original}\t{others}\t{snpeff_ann}')  # snpeff_ann already contains a \n
                    yield [start_original, others, snpeff_ann]

    def call_nucleotide_variants(self, aligner: Callable) -> Generator[Tuple, None, None]:
        """
        Return a generator that iterates over the variants produced by snpEff and each time it yields a tuple of values
        describing a row of the table NucleotideVariant. The tuple reports in order:
        sequence original, sequence alternative, start original, start alternative, variant length, variant type.
        """
        for variant in self._create_or_read_aligned_variants(aligner):
            start_original, others, snpeff_ann = variant
            variant_type, start_alternative, variant_length, sequence_original, sequence_alternative = others.split(',')
            variant_impacts = set()
            for ann in snpeff_ann.split(","):
                s = ann.split("|")
                variant_impacts.add((s[1], s[2], s[3]))
            # assert len(set(variant_impacts)) == len(variant_impacts), 'variant impacts are not unique in each nucleotide variant'
            yield sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type, variant_impacts

    def annotations_and_amino_acid_variants(self, reference_virus_sample) -> [Tuple]:
        """
        :return: a list of annotations for the Annotation and AminoacidVariant table. Each element of the list is a
        tuple of values expressed in this order:
        start, stop, feature_type, gene_name, product, external_reference, amino_acid_sequence, amino_acid_variants
        where amino_acid_variants is either None (if the current sample is the reference sequence for this virus) or
        a list of tuples describing the amino acid variants against the reference sample and described as:
        start position, reference amino acid(s), alternative amino acid(s), variant length, variant type
        """
        if not self._annotations:
            def merge_intervals(e):
                intervals = e.xpath(".//INSDInterval")
                intervals2 = []
                for i in intervals:
                    start = int(text_at_node(i, './/INSDInterval_from'))
                    stop = int(text_at_node(i, './/INSDInterval_to'))
                    intervals2.append((start, stop))

                if intervals:
                    min_start = min(x[0] for x in intervals2)
                    max_stop = max(x[1] for x in intervals2)
                    return min_start, max_stop
                else:
                    return None, None

            def get_annotation(e):
                start, stop = merge_intervals(e)
                feature_type = text_at_node(e, './/INSDFeature_key')
                gene_name = text_at_node(e, './/INSDQualifier[./INSDQualifier_name/text() = "gene"]/INSDQualifier_value',
                                         False)

                product = text_at_node(e, './/INSDQualifier[./INSDQualifier_name/text() = "product"]/INSDQualifier_value',
                                       False)
                db_xref = text_at_node(e, './/INSDQualifier[./INSDQualifier_name/text() = "db_xref"]/INSDQualifier_value',
                                       False)
                protein_id = text_at_node(e,
                                          './/INSDQualifier[./INSDQualifier_name/text() = "protein_id"]/INSDQualifier_value',
                                          False)
                amino_acid_sequence = text_at_node(e,
                                                  './/INSDQualifier[./INSDQualifier_name/text() = "translation"]/INSDQualifier_value',
                                                  False)

                if protein_id:
                    protein_id = "ProteinID:" + protein_id

                # merge with comma
                db_xref_merged = [x for x in [protein_id, db_xref] if x is not None]

                db_xref_merged = ','.join(db_xref_merged)
                #  select one of them:
                #         db_xref_merged = coalesce(db_xref_merged,'db_xref', mandatory=False, multiple=True)

                return start, stop, feature_type, gene_name, product, db_xref_merged, amino_acid_sequence

            annotations_res = []
            for a_feature in self.sample_xml.xpath(".//INSDFeature"):
                try:
                    annotation = get_annotation(a_feature)
                except AssertionError:
                    pass
                else:
                    if annotation and annotation[2] != 'source':
                        aa_variants = self._call_amino_acid_variants(annotation, reference_virus_sample)
                        annotations_res.append((*annotation, aa_variants))
            self._annotations = annotations_res
        return self._annotations

    def _call_amino_acid_variants(self, annotation: Tuple, reference_virus_sample) -> List[Tuple]:
        """
        :param annotation: one of the annotations from VirusSample.annotation()
        :param reference_virus_sample: the VirusSample instance of the reference sample
        :return: given an annotation, return a list of amino acid variants, each one described as a tuple:
        start position, reference amino acid(s), alternative amino acid(s), variant length, variant type
        """
        if reference_virus_sample is None or self == reference_virus_sample:
            return None
        else:
            def call_variants(aa_ref, aa_seq) -> List[Tuple]:
                """
                For every pair of amino acid sequences, returns in order:
                start position, reference amino acid(s), altenative amino acid(s), variant length, variant type
                """
                alignment_aa = pairwise2.align.globalms(aa_ref, aa_seq, 2, -1, -1, -.5)
                ref_aligned_aa = alignment_aa[0][0]
                seq_aligned_aa = alignment_aa[0][1]

                ref_positions_aa = [0 for i in range(len(seq_aligned_aa))]
                pos = 0
                for i in range(len(ref_aligned_aa)):
                    if ref_aligned_aa[i] != '-':
                        pos += 1
                    ref_positions_aa[i] = pos

                aligned_aa = list(zip(ref_positions_aa, ref_aligned_aa, seq_aligned_aa))
                mut_set = set()
                for t in aligned_aa:
                    if t[2] == "-":
                        mutpos = t[0]
                        original = t[1]
                        alternative = t[2]
                        mut_len = max(len(original), len(alternative))
                        mut_type = "DEL"
                        mut_set.add((mutpos, original, alternative, mut_len, mut_type))
                    elif t[1] != t[2]:
                        mutpos = t[0]
                        original = [r for (p, r, s) in aligned_aa if p == mutpos if r != "-"][0]
                        alternative = "".join([s for (p, r, s) in aligned_aa if p == mutpos])
                        mut_len = max(len(original), len(alternative))
                        mut_type = "SUB"
                        mut_set.add((mutpos, original, alternative, mut_len, mut_type))
                    elif t[1] == "-":
                        mutpos = t[0]
                        original = t[1]
                        alternative = t[2]
                        mut_len = max(len(original), len(alternative))
                        mut_type = "SUB"
                        mut_set.add((mutpos, original, alternative, mut_len, mut_type))
                aa_variants = []
                for mut in mut_set:
                    aa_variants.append(mut)

                return aa_variants

            reference_annotations: List[Tuple] = reference_virus_sample.annotations_and_amino_acid_variants(None)
            common_gene_annotations = []
            _, _, _, gene_name, _, _, amino_acid_sequence = annotation
            if gene_name and amino_acid_sequence:
                for _, _, _, ref_gene_name, _, _, ref_amino_acid_sequence, _ in reference_annotations:
                    if ref_gene_name and ref_amino_acid_sequence and ref_gene_name == gene_name:
                        common_gene_annotations.append((ref_amino_acid_sequence, amino_acid_sequence))

            variants = []
            for ref, non_ref in common_gene_annotations:
                for variant in call_variants(ref, non_ref):
                    variants.append(variant)
            return variants



    # def amino_acid_variants(self, reference: List[Tuple[str, Optional[str]]]) -> Generator[Tuple, None, None]:
    #     common_gene_aa = []
    #     for _, _, _, gene_name, _, _, amino_acid_sequence in self.annotations():
    #         for reference_gene_name, reference_aa_sequence in reference:
    #             if gene_name == reference_gene_name:
    #                 common_gene_aa.append((gene_name, amino_acid_sequence, reference_aa_sequence))
    #
    #     def call_aa_variants(aa_ref, aa_seq, name):
    #         aa_variants = []
    #         alignment_aa = pairwise2.align.globalms(aa_ref, aa_seq, 2, -1, -1, -.5)
    #         ref_aligned_aa = alignment_aa[0][0]
    #         seq_aligned_aa = alignment_aa[0][1]
    #
    #         ref_positions_aa = np.zeros(len(seq_aligned_aa), dtype=int)
    #         pos = 0
    #         for i in range(len(ref_aligned_aa)):
    #             if ref_aligned_aa[i] != '-':
    #                 pos += 1
    #             ref_positions_aa[i] = pos
    #
    #         aligned_aa = list(zip(ref_positions_aa, ref_aligned_aa, seq_aligned_aa))
    #         mut_set = set()
    #         for t in aligned_aa:
    #             if t[2] == "-":
    #                 mutpos = t[0]
    #                 original = t[1]
    #                 alternative = t[2]
    #                 mut_type = "DEL"
    #                 mut_set.add((name, mutpos, original, alternative, mut_type))
    #             elif t[1] != t[2]:
    #                 mutpos = t[0]
    #                 original = [r for (p, r, s) in aligned_aa if p == mutpos if r != "-"][0]
    #                 alternative = "".join([s for (p, r, s) in aligned_aa if p == mutpos])
    #                 mut_type = "SUB"
    #                 mut_set.add((name, mutpos, original, alternative, mut_type))
    #             elif t[1] == "-":
    #                 mutpos = t[0]
    #                 original = t[1]
    #                 alternative = t[2]
    #                 mut_type = "SUB"
    #                 mut_set.add((name, mutpos, original, alternative, mut_type))
    #         for mut in mut_set:
    #             aa_variants.append(mut)
    #
    #         return aa_variants
    #
    #     for name, aa_seq, ref_aa_seq in common_gene_aa:
    #         yield call_aa_variants(aa_ref=ref_aa_seq, aa_seq=aa_seq, name=name)

    # def amino_acid_variants(self, reference_virus_sample) -> Generator[Tuple, None, None]:
    #     common_gene_aa = []
    #     for _, _, _, gene_name, _, _, amino_acid_sequence in self.annotations():
    #         for _, _, _, reference_gene_name, _, _, reference_aa_sequence in reference_virus_sample.annotations():
    #             if reference_gene_name and reference_aa_sequence and gene_name == reference_gene_name:
    #                 common_gene_aa.append((gene_name, amino_acid_sequence, reference_aa_sequence))
    #
    #     def call_aa_variants(aa_ref, aa_seq):
    #         aa_variants = []
    #         alignment_aa = pairwise2.align.globalms(aa_ref, aa_seq, 2, -1, -1, -.5)
    #         ref_aligned_aa = alignment_aa[0][0]
    #         seq_aligned_aa = alignment_aa[0][1]
    #
    #         ref_positions_aa = np.zeros(len(seq_aligned_aa), dtype=int)
    #         pos = 0
    #         for i in range(len(ref_aligned_aa)):
    #             if ref_aligned_aa[i] != '-':
    #                 pos += 1
    #             ref_positions_aa[i] = pos
    #
    #         aligned_aa = list(zip(ref_positions_aa, ref_aligned_aa, seq_aligned_aa))
    #         mut_set = set()
    #         for t in aligned_aa:
    #             if t[2] == "-":
    #                 mutpos = t[0]
    #                 original = t[1]
    #                 alternative = t[2]
    #                 mut_type = "DEL"
    #                 mut_set.add((mutpos, original, alternative, mut_type))
    #             elif t[1] != t[2]:
    #                 mutpos = t[0]
    #                 original = [r for (p, r, s) in aligned_aa if p == mutpos if r != "-"][0]
    #                 alternative = "".join([s for (p, r, s) in aligned_aa if p == mutpos])
    #                 mut_type = "SUB"
    #                 mut_set.add((mutpos, original, alternative, mut_type))
    #             elif t[1] == "-":
    #                 mutpos = t[0]
    #                 original = t[1]
    #                 alternative = t[2]
    #                 mut_type = "SUB"
    #                 mut_set.add((mutpos, original, alternative, mut_type))
    #         for mut in mut_set:
    #             aa_variants.append(mut)
    #
    #         return aa_variants
    #
    #     for name, aa_seq, ref_aa_seq in common_gene_aa:
    #         yield call_aa_variants(aa_ref=ref_aa_seq, aa_seq=aa_seq, name=name)


def structured_comment(el, key):
    comment = text_at_node(el, './/INSDSeq_comment' , False)
    if comment:
        sub = re.findall("##Assembly-Data-START## ;(.*); ##Assembly-Data-END##", comment)
        assert len(sub) <=1, f"multiple structured_comment {len(sub)}"
        if sub:
            sub = sub[0]
            subs = sub.split(";")
            subs = [x for x in subs if key in x]
            assert len(subs) <=1, f"multiple structured_comment for key: {key} {len(subs)}"
            if subs:
                return subs[0].split("::")[1].strip()
            else:
                return None
        else:
            return None
    else:
        return None

