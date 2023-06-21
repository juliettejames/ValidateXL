import xml.etree.ElementTree as ET
import pytest

from validatexl.xml_to_df import XMLToDataFrame


@pytest.fixture()
def mz_xml_root():
    """
    Returns test xml file which contains 3 'spectrum_search' attributes, one
    with 3 'search_hit' attributes.
    """
    xml_root = None
    with open("test_data/test.xml", 'r') as xml_file:
        xml_root = ET.XML(xml_file.read())
    return xml_root


def test_loop_through_xml(mz_xml_root):
    """
    Tests that the loop_through_xml returns one search hit when using 
    XML test file.
    """
    xdf = XMLToDataFrame("/data", "my_file.txt")
    actual = xdf._loop_through_xml(mz_xml_root)
    assert len(actual) == 1
