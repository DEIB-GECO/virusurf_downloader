import re
from collections import Counter
from decimal import Decimal

import lxml
from typing import Tuple, Callable, Generator
from Bio import Entrez
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

    def call_variants(self, aligner: Callable) -> Generator[Tuple[str], None, None]:
        """
        Return a generator that iterates over the variants produced by snpEff and each time it yields a tuple of values
        describing a row of the table NucleotideVariant. The tuple reports in order:
        sequence original, sequence alternative, start original, start alternative, variant length, variant type.
        """
        variants: Tuple = aligner(sequence=self.nucleotide_sequence(), sequence_id=self.internal_accession_id)
        snpEff_annotated_variants = variants[0]
        for variant in snpEff_annotated_variants:
            _, start_original, _, _, _, _, others, _ = variant.split("\t")
            variant_type, start_alternative, variant_length, sequence_original, sequence_alternative = others.split(',')
            yield sequence_original, sequence_alternative, start_original, start_alternative, variant_length, variant_type

    def annotations(self) -> [Tuple]:
        """
        :return: a list of annotations for the Annotation table. Each element of the list is a tuple of values expressed
        in this order: start, stop, feature_type, gene_name, product, external_reference, amino_acid_sequence
        """
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
                if annotation:
                    annotations_res.append(annotation)
        return annotations_res


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

