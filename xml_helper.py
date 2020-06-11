from lxml.etree import ElementTree


def text_at_node(root_node: ElementTree, xpath_string, mandatory=True):
    res = root_node.xpath(xpath_string)
    assert len(res) == 1 if mandatory else len(res) <= 1, f"{xpath_string} is available {len(res)} time(s)"

    if res:
        return res[0].text
    else:
        return None


def print_element_tree(etree: ElementTree):
    print(etree.tostring(etree, encoding="unicode", pretty_print=True))

