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

    def _get_xml_root(self):
        """
        Reads the xml_file and returns the root of the xml data using 
        elementtree.

        Returns
        -------
        root: `xml.etree.ElementTree.Element`
            Root node created by ElementTree parser.
        """
        print("Creating dataframe for xml file %s" % self.xml_file)
        xml_data = open(os.path.join(self.xml_path, self.xml_file)).read()
        root = ET.XML(xml_data) # element tree
        return root

    def _get_first_search_hit(self, spec):
        """
        Identifies the first match between a spectra and a crosslink in the 
        xQuest result file. If there is no search_hit an IndexError is raised
        If there is no error raised the first search_hit is returned.

        Parameters
        ----------
        spec: `xml.etree.ElementTree.Element`
            A node representing the spectrum or MS2 scan from the orignal
            Mass spec data file. 

        Returns
        -------
        hit_dict: `dict`
            A dictionary containing the attributes from the search_hit.
        """
        try:
            hit = list(spec.iter('search_hit'))[0]
        except IndexError:
            pass
        else:
            # attributes are stored as dictionary by default
            hit_dict = hit.attrib
            return hit_dict

    def _loop_through_xml(self, root):
        """
        Iterates through the specturm_search object for matches between 
        spectra and crosslinks in the xQuest result file. Calls
        _get_first_search_hit() to store the first search hit for each 
        spectra. Appends the search hit to a dictionary, hit_dict. 
        Appends the hit_dict to a list of all_records.
        
        Parameters
        ----------
        root: `xml.etree.ElementTree.Element`
            Root node created by ElementTree parser.

        Returns
        -------
        all_records: `list`
            List of search_hit records, containing a dictionary of attributes
            associated with each crosslink-spectra match from the xQuest
            result file. 
        """
        all_records = []
        for spec in root.iter('spectrum_search'):
            hit_dict = self._get_first_search_hit(spec)
            if hit_dict:
                all_records.append(hit_dict)
        return all_records

    def _create_dataframe(self, all_records):
        """
        Creates a Pandas DataFrame for all_records. Rows are each crosslink-
        spectra match. Columns are the attributes from the xQuest result file

        Parameters
        ----------
        all_records: `list`
            List of search_hit records, containing a dictionary of attributes
            associated with each crosslink-spectra match from the xQuest
            result file. 

        Returns
        -------
        xml_df: `pd.DataFrame`
            A DataFrame containing crosslink-spectra matches from the xQuest
            result file.
        """
        xml_df = pd.DataFrame(all_records)
        return xml_df

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
        root = self._get_xml_root()
        all_records = self._loop_through_xml(root)
        xml_df = self._create_dataframe(all_records)
        return xml_df
