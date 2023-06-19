import os
import xml.etree.ElementTree as ET

import pandas as pd


class XMLToDataFrame:
    """
    This callable class will generate a Pandas DataFrame which contains
    all the information from the xQuest merged_xml result file.
    A list if DataFrames is returned allowing multiple xQuest result
    files to be parsed.
    
    Parameters
    ----------
    xml_path: `str`
        Path to XML file
    xml_file: `str`
        Name of xQuest XML result file
    """
    def __init__(self, xml_path, xml_file):
        self.xml_path = xml_path
        self.xml_file = xml_file

    def __call__(self):
        """
        Generates a Pandas DataFrame containg the data from the xQuest
        merged_xml result file.

        Returns
        -------
        xml_df: `pd.DataFrame`
            DataFrame containing crosslink identifications from the
            xQuest result file.
        """
        print("Creating dataframe for xml file %s" % self.xml_file)
        xml_data = open(os.path.join(self.xml_path, self.xml_file)).read()
        root = ET.XML(xml_data) # element tree
        
        all_records = []
        for spec in root.iter('spectrum_search'):
            # look through atrrib spec_search for a hit
            # if found take the first (top scoring) hit
            try:
                hit = list(spec.iter('search_hit'))[0]
            except IndexError:
                pass
            else:
                # attributes are stored as dictionary by default
                hit_dict = hit.attrib
                # append dictionary to list
                all_records.append(hit_dict)
        # convert list of dictionaries to dataframe
        xml_df = pd.DataFrame(all_records)
        return xml_df
