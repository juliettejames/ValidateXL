import os

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

from xml_to_df import XMLToDataFrame

class XlSeqId:

    def __init__(self, fasta_path, fasta_name, clean_df):
        self.fasta_file = self._create_fasta_file(fasta_path, fasta_name)
        self.clean_df = clean_df

    def _create_fasta_file(self, fasta_path, fasta_name):
        return os.path.join(fasta_path, fasta_name)
        

    def _extract_amino_acid_position(self, df, record_dict, pep_idx):
        """
        Make lists of the column values for protein id sequence and topology
        loops through searching and extracting amino acid position number using 
        Biopython Seq.find(). Finally adds the position to the topology to get
        AbsPos1 and Abspos2 for both peptides and creates a new DF column.
        This function is called internally by get_seq_id().

        Parameters
        ----------
        df: `pd.DataFrame`
            DataFrame containing the cleaned records from the xQuest xml file
        record_dict: `dict`
            Dictionary of the header details from the FASTA file obtained by 
            Biopython.SeqIO
        pep_idx: `int`
            1 or 2 indicating whether the peptide in the crosslink is the 
            alpha (1) or beta (2) peptide.
        """
        protein = list(df['prot%s' % pep_idx])
        xl = list(df['seq%s' % pep_idx])
        top = list(df['top%s' % pep_idx].apply(pd.to_numeric))
        AbsPos = []
        for i, prot in enumerate(protein):
            my_seq = Seq(str(record_dict[prot].seq))
            pos = my_seq.find(xl[i])
            AbsPos.append(pos + top[i])
        df['AbsPos%s' % pep_idx] = AbsPos


    def get_seq_id(self):
        """
        Makes a dictionary of all proteins in the fasta. Extracts the topology
        of the crosslinker position. Replaces oxidised Methionine in seq1 and 
        seq2 but NOT in Id as this affects the substring matching from the
        fasta file and AbsPos calculation. Calls extract_amino_acid_position
        to find the correct position of the crosslinker.

        Returns
        -------
        abspos_df: `pd.DataFrame`
            DataFrame containing position of the crosslink in the protein
            rather than just the position in the peptide.
        """
        # Make dictionary of all proteins in fastafile, key = xQ protein Id
        # Value = Objects including seq
        record_dict = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))
        print("Obtaining Crosslink Position from Fasta File...")

        # Create Alpha and Beta topology indices
        topolgy = self.clean_df['id'].str.extract(
            "[A-Z]+-[A-Z]+-a(?P<atop>\d+)-b(?P<btop>\d+)", expand=False
        )
        self.clean_df['top1'] = topolgy.atop
        self.clean_df['top2'] = topolgy.btop

        # Replace all X with M in Id as oxidation state of M is irrelevant 
        # for position
        self.clean_df['seq1'] = self.clean_df['seq1'].str.replace('X', 'M')
        self.clean_df['seq2'] = self.clean_df['seq2'].str.replace('X', 'M')

        self._extract_amino_acid_position(self.clean_df, record_dict, 1)
        self._extract_amino_acid_position(self.clean_df, record_dict, 2)

        abspos_df = self.clean_df.copy()
        return abspos_df 

    def __call__(self):
        """
        Generates DataFrame containing xQuest results and the postion of the 
        crosslink in the overall protein sequence.

        Returns
        -------
        abspos_df: `pd.DataFrame`
            DataFrame containing position of the crosslink in the protein
            rather than just the position in the peptide.
        """
        abspos_df = self.get_seq_id()
        return abspos_df
