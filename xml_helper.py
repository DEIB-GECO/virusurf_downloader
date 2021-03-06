from typing import Optional

from lxml.etree import ElementTree, tostring


def text_at_node(root_node: ElementTree, xpath_string, mandatory=True) -> Optional[str]:
    res = root_node.xpath(xpath_string)
    assert len(res) == 1 if mandatory else len(res) <= 1, f"{xpath_string} is available {len(res)} time(s)"

    if res:
        return res[0].text
    else:
        return None


def print_element_tree(node: ElementTree):
    print(tostring(node, encoding="unicode", pretty_print=True))

