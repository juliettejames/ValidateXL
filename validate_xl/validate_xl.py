import os

from Bio import SeqIO
from Bio.Seq import Seq
import click
import pandas as pd

from xml_to_df import XMLToDataFrame
from xl_seq_id import XlSeqId

def clean_xml_data(xml_df):
    """
    Cleans data scraped from xml file to change '-' to zero in count
    of matched ions. Makes numbers stored as strings numeric. Creates new 
    DataFrame with desired columns and returns only crosslinks.
    """
    print("Cleaning data...")
    # Replace "-" with 0 as not correct in xml attrib.tag 
    xml_df['num_of_matched_ions_beta'] = xml_df[
        'num_of_matched_ions_beta'
    ].str.replace('-','0')
    # Convert strings to numbers
    xml_df[[
        'score',
        'num_of_matched_common_ions_alpha',
        'num_of_matched_common_ions_beta',
        'num_of_matched_xlink_ions_alpha',
        'num_of_matched_xlink_ions_beta',
        'num_of_matched_ions_alpha',
        'num_of_matched_ions_beta'
    ]] = xml_df[[
        'score',
        'num_of_matched_common_ions_alpha',
        'num_of_matched_common_ions_beta',
        'num_of_matched_xlink_ions_alpha',
        'num_of_matched_xlink_ions_beta',
        'num_of_matched_ions_alpha',
        'num_of_matched_ions_beta'
    ]].apply(pd.to_numeric)
    xml_df['poss_comm_matches_alpha'] = xml_df['seq1'].str.len()
    xml_df['poss_comm_matches_beta'] = xml_df['seq2'].str.len()

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
        'poss_comm_matches_alpha',
        'poss_comm_matches_beta'
    ]
    # extract desired columns
    ab_ions = xml_df[columns].copy()
    # remove decoys
    tp_ab = ab_ions[
        ~((ab_ions['prot1'].str.contains("decoy") == True) | 
        (ab_ions['prot2'].str.contains("decoy") == True))
    ]
    # extract crosslinks only
    clean_df = tp_ab[tp_ab['type'] == 'xlink']
    return clean_df


def validated_results(abspos_df, path, xml_file):
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
    abspos_df['xl_ion_matches'] = abspos_df.num_of_matched_xlink_ions_alpha + \
        abspos_df.num_of_matched_xlink_ions_beta
    abspos_df['Seq_coverage_alpha'] = abspos_df.num_of_matched_common_ions_alpha / \
        abspos_df.poss_comm_matches_alpha
    abspos_df['Seq_coverage_beta'] = abspos_df.num_of_matched_common_ions_beta / \
        abspos_df.poss_comm_matches_beta
    # groups all matching ids (sensitve to linker position)
    # returns the  highest scoring crosslink for each group
    maxes = abspos_df.groupby("id").score.transform(max)
    top_xl = abspos_df[abspos_df.score == maxes]
    top_xl.to_csv(
        os.path.join(
            path, "%s_full.csv"
            ) % xml_file[:-4], float_format='%.2f'
    )
    validated = top_xl.loc[
        (top_xl['xl_ion_matches'] >= 2) &
        (top_xl['Seq_coverage_alpha'] >= 0.3) &
        (top_xl['Seq_coverage_beta'] >= 0.3)
    ]

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

    rejected = top_xl.loc[
        (top_xl['Seq_coverage_alpha'] < 0.3) &
        (top_xl['Seq_coverage_beta'] < 0.3) &
        (top_xl['xl_ion_matches'] / \
            (top_xl['poss_comm_matches_alpha'] + \
                top_xl['poss_comm_matches_beta']) < 0.3)
    ]

    return validated, manual, rejected


def csv_output(validated_df, manual_df, rejected_df, path, xml_file):
    """
    Generates 3 CSV files for the xQuest XML result file;
    Validate, Manual, Rejected. Column headers are labelled to matching
    the xQuest result file to allow further analysis and comparison.
    """
    header = [
        'Id', 'type', 'prot1', 'prot2', 'Spectrum', 'AbsPos1', 'AbsPos2',
        'measured_mass', 'charge', 'MatchOdds', 'Xcorrx', 'Xcorrb', 
        'WTIC', 'Intsum', 'ld-Score', 'xl_ion_matches',
        'Seq_coverage_alpha', 'Seq_coverage_beta'
    ]
    for df, res_df_type in zip(
        (validated_df, manual_df, rejected_df),
        ("validated", "manual", "rejected")
    ):
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
                res_df_type, xml_file
            )
        )
        df.to_csv(
            os.path.join(
                path, "%s_%s_results.csv"% (
                    xml_file[:-4], res_df_type
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
    xml_file = "9mix_imcs_merged_xquest.xml"

    create_df = XMLToDataFrame(xml_path, xml_file)()
    format_df = clean_xml_data(create_df)

    abspos_df = XlSeqId(fasta_path, fasta_file, format_df)() 
    validated_df, manual_df, rejected_df = validated_results(
        abspos_df, output_dir, xml_file
    )
    csv_output(validated_df, manual_df, rejected_df, output_dir, xml_file)


if __name__ == "__main__":
    cli()
