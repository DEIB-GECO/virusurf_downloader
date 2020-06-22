import re
from collections import Counter
from decimal import Decimal

import lxml
from typing import Tuple
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

