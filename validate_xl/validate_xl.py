import os
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.Seq import Seq
import click
import pandas as pd


def xml2df(xml_path, xmls, energies):
    """
    Creates a list of dataframes from the xQuest result xml files.
    
    Parameters
    ----------
    xml_path: `str`
        Path to XML files
    xmls: `str`
        Name of xQuest XML result files
    energies: `str`
        The energy ramp used in the mass spec analysis

    Returns
    -------
    xml_dfs: `list`
        List of DataFrames containing crosslink identifications from the
        xQuest result files.
    """
    xml_dfs = []
    for i, xml in enumerate(xmls):
        print("Creating dataframe for xml file %s" % energies[i])
        xml_data = open(os.path.join(xml_path, xml)).read()
        root = ET.XML(xml_data) # element tree
        
        all_records = []
        for spec in root.iter('spectrum_search'):
            # look through atrrib spec_search for a hit
            #if found take the first (top scoring) hit
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
        df = pd.DataFrame(all_records)
        # append each dataframe to list of dataframes
        xml_dfs.append(df)

    return xml_dfs


def clean_xml_data(xml_dfs, energies):
    """
    Cleans data scraped from xml file to remove change '-' to zero in count
    of matched ions. Makes numbers stored as strings numeric. Creates new 
    DataFrame with desired columns and returns only crosslinks.
    """
    clean_dfs = []
    for i, df in enumerate(xml_dfs):
        print("Cleaning data for %s" % energies[i])
        # Replace "-" with 0 as not correct in xml attrib.tag 
        df['num_of_matched_ions_beta'] = df[
            'num_of_matched_ions_beta'
        ].str.replace('-','0')
        # Convert strings to numbers
        df[[
            'score',
            'num_of_matched_common_ions_alpha',
            'num_of_matched_common_ions_beta',
            'num_of_matched_xlink_ions_alpha',
            'num_of_matched_xlink_ions_beta',
            'num_of_matched_ions_alpha',
            'num_of_matched_ions_beta'
        ]] = df[[
            'score',
            'num_of_matched_common_ions_alpha',
            'num_of_matched_common_ions_beta',
            'num_of_matched_xlink_ions_alpha',
            'num_of_matched_xlink_ions_beta',
            'num_of_matched_ions_alpha',
            'num_of_matched_ions_beta'
        ]].apply(pd.to_numeric)
        #change this to match formula
        df['poss_comm_matches_alpha'] = df['seq1'].str.len()
        df['poss_comm_matches_beta'] = df['seq2'].str.len()

        columns = [
            'id', 'seq1', 'seq2', 'type', 'prot1', 'prot2',
            'match_odds', 'xcorrb', 'xcorrx', 'wTIC', 'intsum',
            'annotated_spec', 'charge', 'measured_mass', 'score',
            'num_of_matched_common_ions_alpha',
            'num_of_matched_common_ions_beta',
            'num_of_matched_xlink_ions_alpha',
            'num_of_matched_xlink_ions_beta',
            'num_of_matched_ions_alpha',
            'num_of_matched_ions_beta',
            #change this name to remove common
            'poss_comm_matches_alpha',
            'poss_comm_matches_beta'
        ]
        # extract desired columns
        ab_ions = df[columns].copy()
        # remove decoys
        tp_ab = ab_ions[
            ~((ab_ions['prot1'].str.contains("decoy") == True) | 
            (ab_ions['prot2'].str.contains("decoy") == True))
        ]
        # extract crosslinks only
        tp_ab_xl = tp_ab[tp_ab['type'] == 'xlink']
        clean_dfs.append(tp_ab_xl)
    return clean_dfs


def extract_amino_acid_position(df, record_dict, pep_idx):
    """
    Make lists of the column values for protein id sequence and topology
    loops through searching and extracting amino acid position number using 
    Biopython Seq.find(). Finally adds the position to the topology to get
    AbsPos1 and Abspos2 for both peptides and creates a new DF column.
    This function is called internally by get_seq_id().
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


def get_seq_id(clean_dfs, fasta_file):
    """
    Makes a dictionary of all proteins in the fasta. Extracts the topology
    of the crosslinker position. Replaces oxidised Methionine in seq1 and 
    seq2 but NOT in Id as this affects the substring matching from the
    fasta file and AbsPos calculation. Calls extract_amino_acid_position
    to find the correct position of the crosslinker.
    """
    # Make dictionary of all proteins in fastafile, key = xQ protein Id
    # Value = Objects including seq
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    print("Obtaining Crosslink Position from Fasta File")

    # Create Alpha and Beta topology indices
    for df in clean_dfs:
        topolgy = df['id'].str.extract(
            "[A-Z]+-[A-Z]+-a(?P<atop>\d+)-b(?P<btop>\d+)", expand=False
        )
        df['top1'] = topolgy.atop
        df['top2'] = topolgy.btop

        # Replace all X with M in Id as oxidation state of M is irrelevant 
        # for position
        df['seq1'] = df['seq1'].str.replace('X', 'M')
        df['seq2'] = df['seq2'].str.replace('X', 'M')

        extract_amino_acid_position(df, record_dict, 1)
        extract_amino_acid_position(df, record_dict, 2)
    return clean_dfs


def validated_results(df_list, path, energies):
    """
    Validates only the highest scoring
    identification for a particual crosslink based on sequence id. 
    Generates Validated, Manual and Rejected CSVs. 
    Validated XLs have at least 2 crosslinked peaks and 30% sequence
    coverage for fragment ions on BOTH alpha and beta peptides.
    Rejected crosslinks have less that 30% coverage for both alpha, 
    beta and crosslinked ions.
    Crosslinks recommended for Manual validation are the remaining set. 
    i.e those which have >30% sequence coverage on alpha OR beta OR
    at least 30% coverage of the crosslinker ions.
    """
    validated_dfs = []
    manual_dfs = []
    rejected_dfs = []
    for i, df in enumerate(df_list):
        df['xl_ion_matches'] = df.num_of_matched_xlink_ions_alpha + \
            df.num_of_matched_xlink_ions_beta
        df['Seq_coverage_alpha'] = df.num_of_matched_common_ions_alpha / \
            df.poss_comm_matches_alpha
        df['Seq_coverage_beta'] = df.num_of_matched_common_ions_beta / \
            df.poss_comm_matches_beta
        # groups all matching ids (sensitve to linker position)
        # returns hte  higest scoring crosslink for each group
        maxes = df.groupby("id").score.transform(max)
        top_xl = df[df.score == maxes]
        top_xl.to_csv(
            os.path.join(
                path, "%s_full.csv"
            ) % energies[i], float_format='%.2f'
        )
        validated = top_xl.loc[
            (top_xl['xl_ion_matches'] >= 2) &
            (top_xl['Seq_coverage_alpha'] >= 0.3) &
            (top_xl['Seq_coverage_beta'] >= 0.3)
        ]
        validated_dfs.append(validated)

        manual = top_xl.loc[
            (
                (top_xl['Seq_coverage_alpha'] >= 0.3) &
                (top_xl['Seq_coverage_beta'] < 0.3)
            ) | (
                (top_xl['Seq_coverage_alpha'] < 0.3) &
                (top_xl['Seq_coverage_beta'] >= 0.3)
            ) | (
                (top_xl['Seq_coverage_alpha'] < 0.3) &
                (top_xl['Seq_coverage_beta'] < 0.3) &
                (top_xl['xl_ion_matches'] / \
                (top_xl['poss_comm_matches_alpha'] + \
                    top_xl['poss_comm_matches_beta']) >= 0.3))
        ]
        manual_dfs.append(manual)

        rejected = top_xl.loc[
            (top_xl['Seq_coverage_alpha'] < 0.3) &
            (top_xl['Seq_coverage_beta'] < 0.3) &
            (top_xl['xl_ion_matches'] / \
                (top_xl['poss_comm_matches_alpha'] + \
                    top_xl['poss_comm_matches_beta']) < 0.3)
        ]
        rejected_dfs.append(rejected)

    return validated_dfs, manual_dfs, rejected_dfs


def csv_output(validated_dfs, manual_dfs, rejected_dfs, path, energies):
    """
    Generates 3 CSV files for each xQuest XML result file;
    Validate, Manual, Rejected. Column headers are labelled to matching
    the xQuest result file to allow further analysis and comparison.
    """
    header = [
        'Id', 'type', 'prot1', 'prot2', 'Spectrum', 'AbsPos1', 'AbsPos2',
        'measured_mass', 'charge', 'MatchOdds', 'Xcorrx', 'Xcorrb', 
        'WTIC', 'Intsum', 'ld-Score', 'xl_ion_matches',
        'Seq_coverage_alpha', 'Seq_coverage_beta'
    ]
    for dfs, res_df_type in zip(
        (validated_dfs, manual_dfs, rejected_dfs),
        ("validated", "manual", "rejected")
    ):
        for i, df in enumerate(dfs):
            df.rename(
                columns={
                    'score': 'ld-Score',
                    'id': 'Id',
                    'wTIC': 'WTIC',
                    'xcorrb': 'Xcorrb',
                    'xcorrx': 'Xcorrx',
                    'match_odds': 'MatchOdds',
                    'intsum': 'Intsum',
                    'annotated_spec': 'Spectrum',
                }, inplace=True
            )
            print(
                "Generating %s results for %s" % (
                    res_df_type, energies[i]
                )
            )
            df.to_csv(
                os.path.join(
                    path, "%s_%s_results.csv"% (
                        energies[i], res_df_type
                    )
                ), columns=header, float_format='%.2f'
            )


@click.command()
@click.option('--input-dir', 'input_dir', default=None, help='The location of the Fasta and xQuest merged_xml input files')
@click.option('--output-dir', 'output_dir', default=None, help='The location to output the CSV result files')
def cli(input_dir, output_dir):
    fasta_path = os.path.join(input_dir, "fasta")
    fasta_file = os.path.join(fasta_path, "9_mix.fasta")
    xml_path = os.path.join(input_dir, "xq_xmls")

    energies = ['mid']
    xmls = [
        "9mix_imcs_merged_xquest.xml" 
    ]

    create_df = xml2df(xml_path, xmls, energies)
    format_dfs = clean_xml_data(create_df, energies)
    abspos_dfs = get_seq_id(format_dfs, fasta_file)
    validated_dfs, manual_dfs, rejected_dfs = validated_results(
        abspos_dfs, output_dir, energies
    )
    csv_output(validated_dfs, manual_dfs, rejected_dfs, output_dir, energies)


if __name__ == "__main__":
    cli()
