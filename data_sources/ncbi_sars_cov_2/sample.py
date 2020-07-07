import re
from collections import Counter
from decimal import Decimal
from lxml.etree import ElementTree

from typing import Tuple, Callable, Generator, List, Iterable, Optional
from Bio import pairwise2, Entrez
from loguru import logger
from lxml import etree

import IlCodiCE
from geo_groups import geo_groups
from locations import *
from data_sources.virus_sample import VirusSample
from xml_helper import text_at_node
from database_tom import RollbackTransactionWithoutError
from datetime import datetime
from dateutil.parser import parse


default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)


class NCBISarsCov2Sample(VirusSample):
    virus_name = 'NCBI_sars_cov_2'
    cached_taxon_id = {}
    re_anna_RuL3z = re.compile(r'(\d+)([,|.]?)(\d*).*')

    """
    Set of getters to take the relevant information from a virus sample in INSDSeq XML format.
    Example of virus sample @ https://www.ncbi.nlm.nih.gov/nuccore/MN908947
    INSDSeq XML format DTD @ https://www.ncbi.nlm.nih.gov/data_specs/dtd/INSD_INSDSeq.mod.dtd
    """
    def __init__(self, virus_sample_file_path: str, internal_accession_id):
        """
        :param virus_sample_file_path: virus sample file (INDSeq XML) path
        :param internal_accession_id: sequence accession id ( == numeric part of a GI id, e.g. '1798174254' of 'gi|1798174254')
        Raises lxml.etree.XMLSyntaxError if the source file is malformed or empty.
        """
        super().__init__()
        self.sample_xml: ElementTree = etree.parse(virus_sample_file_path, parser=etree.XMLParser(remove_blank_text=True))
        self.internal_accession_id = internal_accession_id

        # init internal variables
        self._journal = None
        self._host_value_list = None
        self._annotations = None
        self._nuc_seq = None

    def internal_id(self):
        return self.internal_accession_id

    def primary_accession_number(self):
        return text_at_node(self.sample_xml, './/INSDSeq_accession-version')

    def alternative_accession_number(self):
        return str(self.internal_accession_id) if self.internal_accession_id else None

    def strain(self):
        strain = text_at_node(self.sample_xml, './/INSDQualifier[./INSDQualifier_name/text() = "strain"]/INSDQualifier_value', False)
        if not strain:
            # try to get the strain from the filed FEATURES -> SOURCE -> /isolate
            strain = text_at_node(self.sample_xml, './/INSDQualifier[./INSDQualifier_name/text() = "isolate"]/INSDQualifier_value', False)
        return strain

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
            length = self.length()
            if length and length < 28407:      # 95 % of the length of the reference sequence
                return False
            else:
                return None

    def nucleotide_sequence(self):
        if not self._nuc_seq:
            try:
                self._nuc_seq = text_at_node(self.sample_xml, './/INSDSeq_sequence')
            except AssertionError:
                raise RollbackTransactionWithoutError(f'Insertion of sequence {self.internal_accession_id} skipped because of missing nucleotide seq.')
        return self._nuc_seq

    # noinspection PyMethodMayBeStatic
    def strand(self):
        return 'positive'

    def length(self):
        length = int(text_at_node(self.sample_xml, './/INSDSeq_length'))
        assert length == len(self.nucleotide_sequence())
        return length

    def gc_percent(self):
        c = Counter(self.nucleotide_sequence().lower())
        count_known_nucleotides = (c['g'] + c['c'] + c['a'] + c['t'] + c['u'])
        if count_known_nucleotides != 0:
            gc_percentage = (c['g'] + c['c']) / count_known_nucleotides * 100
            gc_percentage = Decimal(gc_percentage)
            gc_percentage = round(gc_percentage, 2)
            return gc_percentage
        else:
            return 0

    def n_percent(self):
        length = self.length()
        if length != 0:
            c = Counter(self.nucleotide_sequence().lower())
            n_percentage = (c['n']) / length * 100
            n_percentage = Decimal(n_percentage)
            n_percentage = round(n_percentage, 2)
            return n_percentage
        else:
            return 0

    def lineage(self):
        return None

    def clade(self):
        return None

    def sequencing_technology(self):
        return _structured_comment(self.sample_xml, 'Sequencing Technology')

    def assembly_method(self):
        return _structured_comment(self.sample_xml, 'Assembly Method')

    def coverage(self):
        _input = _structured_comment(self.sample_xml, 'Coverage')
        if not _input:
            return None
        info = self.re_anna_RuL3z.match(_input)
        if not info:
            raise ValueError(f'coverage string {_input} doesn\'t match the regex.')
        else:
            o1 = info.group(1)  # integer part
            o3 = info.group(3)  # decimal (optional)
            # apply Anna's rules
            output = int(o1)
            if o3 and len(o3) < 3:
                decimals = round(int(o3) / 10)
                output += decimals
            return output

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
        node = text_at_node(self.sample_xml,
                            '..//INSDQualifier[./INSDQualifier_name/text() = "country"]/INSDQualifier_value',
                            mandatory=False)
        if node:
            node = node.split(":")
            country = node[0]
            region = node[1] if len(node) > 1 else None
        geo_group = geo_groups.get(country.lower()) if country else None
        return country, region, geo_group

    def _init_and_get_journal(self):
        if not self._journal:
            references = self.sample_xml.xpath('.//INSDReference[./INSDReference_title/text() = "Direct Submission"]')
            assert len(references) > 0, 'there must be at least one direct submission'
            reference = references[0]
            journal_string = text_at_node(reference, "./INSDReference_journal")
            assert journal_string.startswith(
                "Submitted "), 'Cannot find submitted in the Journal of direct submission reference'
            self._journal = re.split("[()]", journal_string, maxsplit=2)[1:]
            assert len(self._journal) == 2, f"Journal value problem '{journal_string}' {self._journal}"

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

    def submission_date(self):
        try:
            return datetime.strptime(self._init_and_get_journal()[0], '%d-%b-%Y')
        except TypeError:
            return None

    def originating_lab(self) -> Optional[str]:
        return None

    def sequencing_lab(self) -> Optional[str]:
        try:
            return self._init_and_get_journal()[1]
        except TypeError:
            return None

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
        taxon_name = self.taxon_name()
        if not taxon_name:
            return None
        else:
            taxon_id = NCBISarsCov2Sample.cached_taxon_id.get(taxon_name)
            if taxon_id == -1:  # -1 means the cached taxon_id for this taxon name was searched before
                return None
            elif taxon_id is None:
                try:
                    with Entrez.esearch(db="taxonomy", term=f'"{taxon_name}"[Scientific Name]', rettype=None,
                                        retmode="xml") as handle:
                        response = Entrez.read(handle)
                        if response['Count'] == '1':
                            taxon_id = int(response['IdList'][0])
                            NCBISarsCov2Sample.cached_taxon_id[taxon_name] = taxon_id
                        else:
                            logger.warning(f'can\'t find the taxon id for taxon name {taxon_name}')
                            NCBISarsCov2Sample.cached_taxon_id[taxon_name] = -1  # save -1 in cache to distinguish from non cached taxon_ids
                            taxon_id = None
                except:
                    logger.exception(f'Exception occurred while fetching the taxon id of {taxon_name} in sample {self.internal_accession_id()}')
            return taxon_id

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

    def _call_or_read_nuc_variants_and_effects(self, aligner: Callable) -> Generator[Iterable, None, None]:
        cache_file_path = f'{get_local_folder_for(source_name=self.virus_name, _type=FileType.NucleotideVariants)}{os.path.sep}{self.internal_accession_id}'
        try:
            with open(cache_file_path, mode='r') as f:
                logger.trace(f'reading nucleotide variants for sample {self.internal_accession_id} from disk')
                for line in f:
                    yield line.rstrip().split('\t')
        except FileNotFoundError:
            # compute them and write to file
            logger.trace(f'calling nucleotide variants for sample {self.internal_accession_id}...')
            variants: Tuple = aligner(sequence=self.nucleotide_sequence(), sequence_id=self.internal_accession_id)
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
        """
        :return: a list of annotations for the Annotation and AminoacidVariant table. Each element of the list is a
        tuple of values expressed in this order:
        start, stop, feature_type, gene_name, product, external_reference, amino_acid_sequence, amino_acid_variants
        where amino_acid_variants is either None (if the current sample is the reference sequence for this virus) or
        a list of tuples describing the amino acid variants against the reference sample and described as:
        start position, reference amino acid(s), alternative amino acid(s), variant length, variant type
        """
        if not self._annotations:
            self._annotations = []
            for a_feature in self.sample_xml.xpath(".//INSDFeature"):
                # compute annotation
                try:
                    start, stop = _merge_intervals(a_feature)
                    feature_type = text_at_node(a_feature, './/INSDFeature_key')
                    if feature_type == 'source':
                        continue
                    gene_name = text_at_node(
                        a_feature,
                        './/INSDQualifier[./INSDQualifier_name/text() = "gene"]/INSDQualifier_value',
                        False)
                    product = text_at_node(
                        a_feature,
                        './/INSDQualifier[./INSDQualifier_name/text() = "product"]/INSDQualifier_value',
                        False)
                    amino_acid_sequence = text_at_node(
                        a_feature,
                        './/INSDQualifier[./INSDQualifier_name/text() = "translation"]/INSDQualifier_value',
                        False)
                    db_xref = text_at_node(
                        a_feature,
                        './/INSDQualifier[./INSDQualifier_name/text() = "db_xref"]/INSDQualifier_value',
                        False)
                    protein_id = text_at_node(
                        a_feature,
                        './/INSDQualifier[./INSDQualifier_name/text() = "protein_id"]/INSDQualifier_value',
                        False)
                    if protein_id:
                        protein_id = "ProteinID:" + protein_id
                    external_reference = [x for x in [protein_id, db_xref] if x is not None]
                    external_reference = ','.join(external_reference) if external_reference else None
                    #  select one of them:
                    #         db_xref_merged = coalesce(db_xref_merged,'db_xref', mandatory=False, multiple=True)
                except AssertionError:
                    continue
                # compute amino acid variants only on CDS except those of the reference
                else:
                    if reference_virus_sample is None or self == reference_virus_sample or not gene_name or not amino_acid_sequence:
                        aa_variants = None
                    else:
                        # ignore annotations named as orf1... and length < 20 k
                        gene_name_lower = gene_name.lower()
                        if gene_name_lower.startswith('orf1') and stop-start < 20000:
                            continue    # skip this annotation

                        # match this CDS with one from the reference
                        common_gene_annotations = []
                        for _, _, _, ref_gene_name_lower, _, _, ref_aa_seq, _ in reference_virus_sample.get_cached_reference_CDS_annotations():
                            if ref_gene_name_lower == gene_name_lower or \
                                    ref_gene_name_lower.startswith('orf1') and gene_name_lower.startswith('orf1'):
                                common_gene_annotations.append((ref_aa_seq, amino_acid_sequence))

                        # call amino acid variants on the pair of amino acid sequences
                        aa_variants = [variant for (aa_ref, aa_seq) in common_gene_annotations for variant in _call_aa_variants(aa_ref, aa_seq)]

                    # append annotation and its amino acid variants
                    annotation_and_aa_variants = (start, stop, feature_type, gene_name, product, external_reference, amino_acid_sequence, aa_variants)
                    self._annotations.append(annotation_and_aa_variants)
                    yield annotation_and_aa_variants
        else:
            for el in self._annotations:
                yield el

    # noinspection PyPep8Naming
    def get_cached_reference_CDS_annotations(self):
        if not hasattr(self, 'cached_CDS_annotations'):
            self.cached_CDS_annotations = []
            for start, stop, feature_type, gene_name, product, external_reference, amino_acid_sequence, aa_variants in self.annotations_and_amino_acid_variants(None):
                if gene_name and amino_acid_sequence:
                    self.cached_CDS_annotations.append((start, stop, feature_type, gene_name.lower(), product, external_reference, amino_acid_sequence, aa_variants))
        return self.cached_CDS_annotations

    def on_before_multiprocessing(self):
        self.nucleotide_sequence()  # make sure nucleotide variants are cached in this sample
        self.sample_xml = None  # release xml reference to prevent errors with multiprocessing

    def nucleotide_var_aligner(self):
        return IlCodiCE.create_aligner_to_reference(reference=self.nucleotide_sequence(),
                                                    annotation_file='sars_cov_2_annotations.tsv',
                                                    is_gisaid=False)


def _structured_comment(el, key):
    comment = text_at_node(el, './/INSDSeq_comment', False)
    if comment:
        sub = re.findall("##Assembly-Data-START## ;(.*); ##Assembly-Data-END##", comment)
        assert len(sub) <= 1, f"multiple structured_comment {len(sub)}"
        if sub:
            sub = sub[0]
            subs = sub.split(";")
            subs = [x for x in subs if key in x]
            assert len(subs) <= 1, f"multiple structured_comment for key: {key} {len(subs)}"
            if subs:
                return subs[0].split("::")[1].strip()
            else:
                return None
        else:
            return None
    else:
        return None


def _merge_intervals(e):
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


def _call_aa_variants(aa_ref, aa_seq) -> List[Tuple]:
    """
    For every pair of amino acid sequences, returns in order:
    start position, reference amino acid(s), altenative amino acid(s), variant length, variant type
    """
    alignment_aa = pairwise2.align.globalms(aa_ref, aa_seq, 2, -1, -1, -.5, one_alignment_only=True)
    ref_aligned_aa = alignment_aa[0][0]
    seq_aligned_aa = alignment_aa[0][1]

    ref_positions_aa = [0 for _ in range(len(seq_aligned_aa))]
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
            mut_set.add((original, alternative, mutpos, mut_len, mut_type))
        elif t[1] != t[2]:
            mutpos = t[0]
            try:
                original = [r for (p, r, s) in aligned_aa if p == mutpos if r != "-"][0]
            except IndexError:
                continue
            alternative = "".join([s for (p, r, s) in aligned_aa if p == mutpos])
            mut_len = max(len(original), len(alternative))
            mut_type = "SUB"
            mut_set.add((original, alternative, mutpos, mut_len, mut_type))
        elif t[1] == "-":
            mutpos = t[0]
            original = t[1]
            alternative = t[2]
            mut_len = max(len(original), len(alternative))
            mut_type = "SUB"
            mut_set.add((original, alternative, mutpos, mut_len, mut_type))
    aa_variants = []
    for mut in mut_set:
        aa_variants.append(mut)

    return aa_variants
