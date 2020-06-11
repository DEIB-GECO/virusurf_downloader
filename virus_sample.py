import re
from collections import Counter
from decimal import Decimal
from Bio import Entrez
from loguru import logger
from lxml import etree
from locations import *
from xml_helper import *


def download_virus_sample_as_xml(sample_accession_id: int) -> str:
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
        self.sample_xml: ElementTree = etree.parse(virus_sample_file_path, parser=etree.XMLParser(remove_blank_text=True))
        self.internal_accession_id = internal_accession_id

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
        return text_at_node(self.sample_xml, './/INSDSeq_sequence')

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

