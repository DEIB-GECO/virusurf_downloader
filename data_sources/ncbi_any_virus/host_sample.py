"""
Created by tomalf2 on nov, 2020.
"""
from datetime import datetime
import lxml
from loguru import logger
from lxml import etree
from typing import Optional, Tuple, List, Collection
from data_sources.ncbi_services import download_or_get_ncbi_host_sample_as_xml
from locations import FileType, get_local_folder_for
from geo_groups import geo_groups


class NCBIHostSample:
    """
    Class dedicated to the parsing of XML files representing entities of the NCBI biosample database.
    If it is possible to obtain these informatin in other ways, avoid using this class as it requires two additional API
    calls for every biosample, thus causing an important performance slowdown measured in an interval ranging from
    +1 second to +6 seconds per sequence, with an average of 4,7 seconds per sequence
    """
    def __init__(self, host_sample_accession_id: str):
        self.acc_id = host_sample_accession_id
        containing_dir = get_local_folder_for('test', FileType.SequenceOrSampleData)
        file_path = download_or_get_ncbi_host_sample_as_xml(containing_dir, self.acc_id)
        self.host_xml: lxml.etree.ElementTree = etree.parse(source=file_path,
                                                              parser=etree.XMLParser(remove_blank_text=True))
        # with open(file_path) as f:
        #     for l in f.readlines():
        #         print(l)

        attribute_nodes_cleaned = []
        attribute_nodes = self.host_xml.xpath('.//Attributes/Attribute')
        # remove attributes with not applicable value
        for node in attribute_nodes:
            value = node.text.lower()
            if 'not' not in value:  # like "not applicable"
                attribute_nodes_cleaned.append(node)
        # save attributes ina key-value pairs
        self.attributes = {}
        for node in attribute_nodes_cleaned:
            if 'attribute_name' in node.attrib:
                key = node.attrib['attribute_name'].lower()
            else:
                raise ValueError('NCBI XML host data is not encoded in any know format')
            value = node.text
            # insert key-value pair or change key
            if key in self.attributes.keys():
                if 'harmonized_name' in node.attrib:
                    key = node.attrib['harmonized_name'].lower()
                    if key in self.attributes.keys():
                        if 'display_name' in node.attrib:
                            key = node.attrib['display_name'].lower()
                            if key in self.attributes.keys():
                                raise ValueError('NCBI XML host data is not encoded in any know format')
            self.attributes[key] = value
        # print(self.attributes)

    def isolation_source(self):
        source = _find_in_attributes_(self.attributes, ('source', 'tissue'))[1]
        if source is not None and 'not' in source:  # like 'not collected'
            source = None
        return source

    def country__region__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        country: Optional[str] = _find_in_attributes(self.attributes, 'country')[1]
        region: Optional[str] = _find_in_attributes(self.attributes, 'region')[1]
        if country is not None:
            if 'not' in country: # like 'not collected'
                country = None
                geo_group = None
            else:
                country = country.strip()
                geo_group = geo_groups.get(country) # up to here country is lowercase
                country = country.capitalize()
        else:
            geo_group = None
        if region is not None:
            if 'not' in region: # like 'not collected'
                region = None
            else:
                region = region.strip().capitalize()
        return country, region, geo_group

    def coverage(self) -> Optional[int]:
        value = _find_in_attributes(self.attributes, 'coverage')[1]
        if value is not None:
            try:
                # noinspection PyTypeChecker
                return round(float(value))
            except ValueError:
                logger.error(f'Error while parsing coverage string {value} from host XML data')
                return None

    def originating_lab(self) -> Optional[str]:
        lab = _find_in_attributes(self.attributes, 'collecting institu')[1]
        if lab is not None and 'not' in lab: # like 'not provided'
            lab = None
        return lab

    def collection_date(self) -> Optional[str]:
        return _find_in_attributes(self.attributes, 'collection date')[1]

    def submission_date(self) -> Optional[datetime]:
        date = _find_in_attributes(self.attributes, 'receipt date')[1]
        if date:
            try:
                return datetime.strptime(date, '%Y-%m-%d')
            except ValueError as e:
                logger.error(e.args)
            except TypeError:
                return None
        else:
            return None

    def host_taxon_name(self) -> Optional[str]:
        host_name = _find_in_attributes(self.attributes, 'scientific')[1]
        if not host_name:
            host_name = _find_in_attributes(self.attributes, 'host')[1]
        return host_name

    def age(self) -> Optional[str]:
        age_value = _find_in_attributes_(self.attributes, ('age', 'years'), ('coverage', 'stage', 'passage'))[1]
        # parse to int to eliminate possible decimals
        if age_value is not None:
            if 'not' in age_value: # like 'not collected'
                age_value = None
            else:
                age_value = _find_str_of_integers(age_value)
        return age_value

    def gender(self) -> Optional[int]:
        gender = _find_in_attributes_(self.attributes, ('sex', 'gender'))[1]
        if gender is not None:
            if 'not' in gender or 'restricted' in gender: # like 'not provided' or 'restricted access'
                gender = None
        return gender


def _find_str_of_integers(string: str) -> Optional[str]:
    string_length = len(string)
    number_found = False
    start_num_idx = 0
    while start_num_idx < string_length and not number_found:
        if string[start_num_idx].isdigit():
            number_found = True
        start_num_idx += 1
    end_num_idx = start_num_idx
    start_num_idx -= 1
    while end_num_idx < string_length and string[end_num_idx].isdigit():
        end_num_idx += 1
    if number_found:
        return string[start_num_idx:end_num_idx]
    else:
        return None


def _find_in_attributes(attributes: dict, keyword, exclude_keyword: Optional[str] = None):
    key = [k for k in attributes if keyword in k]
    if exclude_keyword is not None:
        key = [k for k in key if exclude_keyword not in k]

    if len(key) > 0:
        assert len(key) == 1, f'NCBI XML host data contains more than one valid {[k for k in key]} information'
        value = attributes[key[0]].lower()
        return key[0], value
    else:
        return None, None


def _find_in_attributes_(attributes: dict, include_keywords: Collection[str], exclude_keywords: Optional[Collection[str]] = ()) -> Tuple[Optional[str], Optional[str]]:
    # find eligible attribute keys
    eligible_keys_1 = []
    for word in include_keywords:
        for k in attributes.keys():
            if word in k:
                eligible_keys_1.append(k)

    eligible_keys_2 = []
    # ignore excluded attribute keywords
    for k in eligible_keys_1:
        excluded = False
        for word in exclude_keywords:
            if word in k:
                excluded = True
                break
        if not excluded:
            eligible_keys_2.append(k)

    if len(eligible_keys_2) > 0:
        assert len(eligible_keys_2) == 1, f'NCBI XML host data contains more than one valid {[k for k in eligible_keys_2]} information'
        value = attributes[eligible_keys_2[0]].lower()
        return eligible_keys_2[0], value
    else:
        return None, None
