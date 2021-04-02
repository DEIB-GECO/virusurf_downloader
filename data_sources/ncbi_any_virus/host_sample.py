"""
Created by tomalf2 on nov, 2020.
"""
from datetime import datetime
import lxml
from loguru import logger
from lxml import etree
from typing import Optional, Tuple, Collection
from data_sources.ncbi_services import download_or_get_ncbi_host_sample_as_xml
from locations import remove_file
from geo_groups import geo_groups
import dateutil.parser as dateparser


class NCBIHostSample:
    """
    Class dedicated to the parsing of XML files representing entities of the NCBI biosample database.
    If it is possible to obtain these informatin in other ways, avoid using this class as it requires two additional API
    calls for every biosample, thus causing an important performance slowdown measured in an interval ranging from
    +1 second to +6 seconds per sequence, with an average of 4,7 seconds per sequence
    """

    default_datetime = datetime(2020, 1, 1, 0, 0, 0, 0, None)

    def __init__(self, host_sample_accession_id: str, download_dir: str):
        self.acc_id = host_sample_accession_id
        self.file_path = download_or_get_ncbi_host_sample_as_xml(download_dir, self.acc_id)
        self.host_xml: lxml.etree.ElementTree = etree.parse(source=self.file_path,
                                                            parser=etree.XMLParser(remove_blank_text=True))
        # with open(file_path) as f:
        #     for l in f.readlines():
        #         print(l)

        attribute_nodes_cleaned = []
        attribute_nodes = self.host_xml.xpath('.//Attributes/Attribute')
        # remove attributes with not applicable value
        for node in attribute_nodes:
            value = node.text.lower()
            # ignore values like "not applicable"/"not collected"/"restricted access"/"missing" values
            if 'not' not in value and 'missing' not in value and 'restricted' not in value:
                attribute_nodes_cleaned.append(node)
        # save attributes ina key-value pairs
        self.attributes = {}
        for node in attribute_nodes_cleaned:
            if 'attribute_name' in node.attrib:
                key = node.attrib['attribute_name'].lower()
            else:
                remove_file(self.file_path)
                raise ValueError(f'NCBI XML host {self.acc_id} is not encoded in any know format. XML host file removed.')
            value = node.text
            # insert key-value pair or change key
            if key in self.attributes.keys():
                if 'harmonized_name' in node.attrib:
                    key = node.attrib['harmonized_name'].lower()
                    if key in self.attributes.keys():
                        if 'display_name' in node.attrib:
                            key = node.attrib['display_name'].lower()
                            if key in self.attributes.keys():
                                remove_file(self.file_path)
                                raise ValueError(f'NCBI XML host {self.acc_id} is not encoded in any know format. XML host file removed.')
            self.attributes[key] = value
        # print(self.attributes)

    def isolation_source(self):
        return self._find_in_attributes_(('tissue', 'source'), use_keyword_priority=True)[1]

    def province__region__country__geo_group(self) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        country: Optional[str] = self._find_in_attributes('country')[1]
        region: Optional[str] = self._find_in_attributes('region')[1]
        if country is not None:
            if 'not' in country:  # like 'not collected'
                country = None
                geo_group = None
            else:
                country = country.strip()
                geo_group = geo_groups.get(country)  # up to here country is lowercase
                country = country.capitalize()
        else:
            geo_group = None
        if region is not None:
            if 'not' in region:  # like 'not collected'
                region = None
            else:
                region = region.strip().capitalize()
        return None, region, country, geo_group

    def coverage(self) -> Optional[int]:
        value = self._find_in_attributes('coverage')[1]
        if value is not None:
            try:
                # noinspection PyTypeChecker
                return round(float(value))
            except ValueError:
                remove_file(self.file_path)
                logger.error(f'Error while parsing coverage string {value} from host XML {self.acc_id}. File removed.')
                return None

    def originating_lab(self) -> Optional[str]:
        return self._find_in_attributes('collecting institu')[1]

    def collection_date(self) -> Tuple[Optional[str], Optional[int]]:
        col_d = self._find_in_attributes('collection date')[1]
        if col_d:
            col_d = col_d.strip()
            # find precision (year = 0, month = 1, day = 2
            if '-' in col_d:
                precision = col_d.count('-')
            elif len(col_d) == 4:
                precision = 0
            elif '/' in col_d:
                precision = col_d.count('/')
            elif '\\' in col_d:
                precision = col_d.count('\\')
            elif ' ' in col_d:
                precision = col_d.count(' ')
            else:
                raise AssertionError(
                    f'Unable to parse date string {col_d}. Unexpected separator in date string or no separator.')
            return dateparser.parse(col_d, default=self.default_datetime).strftime('%Y-%m-%d'), precision
        else:
            return None, None

    def submission_date(self) -> Optional[datetime]:
        date = self._find_in_attributes('receipt date')[1]
        if date:
            try:
                return datetime.strptime(date, '%Y-%m-%d')
            except ValueError as e:
                remove_file(self.file_path)
                logger.error(f"XML host file {self.acc_id} removed. {e.args}")
            except TypeError:
                return None
        else:
            return None

    def host_taxon_name(self) -> Optional[str]:
        host_name = self._find_in_attributes('scientific')[1]
        if not host_name:
            host_name = self._find_in_attributes_(('host',), ('disease', 'description'))[1]
        return host_name

    def age(self) -> Optional[str]:
        age_value = self._find_in_attributes_(('age', 'years'), ('coverage', 'stage', 'passage'), True)[1]
        # parse to int to eliminate possible decimals
        if age_value is not None:
            age_value = self._find_str_of_integers(age_value)
        return age_value

    def gender(self) -> Optional[int]:
        return self._find_in_attributes_(('sex', 'gender'))[1]

    @staticmethod
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

    def _find_in_attributes(self, keyword, exclude_keyword: Optional[str] = None):
        key = [k for k in self.attributes if keyword in k]
        if exclude_keyword is not None:
            key = [k for k in key if exclude_keyword not in k]

        if len(key) == 1:
            value = self.attributes[key[0]].lower()
            return key[0], value
        elif len(key) > 1:
            remove_file(self.file_path)
            raise AssertionError(
                f'NCBI XML host {self.acc_id} contains more than one valid {[k for k in key]} information. XML host file removed.')
        else:
            return None, None

    def _find_in_attributes_(self, include_keywords: Collection[str], exclude_keywords: Optional[Collection[str]] = (),
                             use_keyword_priority: bool = False) -> Tuple[Optional[str], Optional[str]]:
        # find eligible attribute keys
        eligible_keys_1 = []
        for word in include_keywords:
            for k in self.attributes.keys():
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

        if len(eligible_keys_2) == 1 or (len(eligible_keys_2) > 1 and use_keyword_priority):
            # eligible_keys are already ordered in the same order as include_keywords
            value = self.attributes[eligible_keys_2[0]].lower()
            return eligible_keys_2[0], value
        elif len(eligible_keys_2) > 1:
            remove_file(self.file_path)
            raise AssertionError(
                f'NCBI XML host {self.acc_id} contains more than one valid {[k for k in eligible_keys_2]} information. XML host '
                f'file removed.')
        else:
            return None, None
