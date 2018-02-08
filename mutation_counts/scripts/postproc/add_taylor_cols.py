#!/usr/bin/env python2

"""
SYNOPSIS

    Annotate the Mutation Hotspots Table with hotspot information from the
    the Chang et al. (2016) paper from Barry Taylor lab (MSK).

NOTES

    (1) Add the following annotation 2 columns:

        (A) "IS_TAYLOR_HOTSPOT": Indicate whether each hotspot was previously
            reported as hotspot the Chang et al. (2016) paper

        (B) "IS_TAYLOR_FALSE_POSITIVE": Indicate whether each hotspot was
            reported as a 'False Positive' by the Chang et al.(2016) paper.

    (2) Chang et al. (2016) Paper is a Codon-level Hotspots Analysis using
        Pan-Cancer TCGA data.

            Paper URL: https://www.ncbi.nlm.nih.gov/pubmed/26619011

    (3) The 'taylor-supplementary_tables_1-5_clean_footer.xls' is a supplemental
        excel spreadsheet provided to us by the Taylor Lab. The footer of the
        sheets was cleaned to allow for easy import.

EXAMPLES

    (1) All Private + All Public
    ./add_taylor_cols.py \
        --in_table ../private/reliance_point_outputs/05May2017/05May2017_ReliancePoint_ClinicalHotspotsTable.txt \
        --out_table tmp/05May2017_ReliancePoint_ClinicalHotspotsTable.CompareTaylor.txt

    (2) Public Only
    ./add_taylor_cols.py \
        --in_table ../private/count_tables/04May2017.count_table.merged.only_public.use_hotspot_stats_model.with_breakdown.txt \
        --out_table tmp/04May2017.count_table.merged.only_public.use_hotspot_stats_model.with_breakdown.CompareTaylor.txt

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import argparse
import pandas as pd

filename = 'taylor-supplementary_tables_1-5_clean_footer.xls'
TAYLOR_SI_TABLES = os.path.abspath(os.path.join(os.path.dirname(__file__), '.', filename))

def main(options):

    # Import Hotspots data from SI Table of Chang at el. (2016) Paper.
    taylor_hotspots_df = pd.read_excel(TAYLOR_SI_TABLES, sheetname="Table S2", skiprows=2)
    taylor_false_positives_df = pd.read_excel(TAYLOR_SI_TABLES, sheetname="Table S3",  skiprows=2)

    # Import the input Hotspots Table
    hotspots_df, in_comment_lines = import_hotspots_table(options.in_table)

    # Set Index of all dfs to be <Gene_Symbol>_<Ref_Codon> (e.g. 'PIK3CA_H1047')
    taylor_hotspots_df.index = taylor_hotspots_df["Hugo Symbol"] + "_" + taylor_hotspots_df["Codon"]

    taylor_false_positives_df.index = taylor_false_positives_df["Hugo Symbol"] +  "_" +  taylor_false_positives_df["Codon"]

    hotspots_df.index = hotspots_df['Gene_Symbol'] + '_' + hotspots_df['Ref_Codon']

    # Consistency Checks
    assert taylor_hotspots_df.index.is_unique
    assert taylor_false_positives_df.index.is_unique
    assert hotspots_df.index.is_unique

    hotspots_df["IS_TAYLOR_HOTSPOT"] = hotspots_df.index.isin(taylor_hotspots_df.index)
    hotspots_df["IS_TAYLOR_FALSE_POSITIVE"] = hotspots_df.index.isin(taylor_false_positives_df.index)

    # Parse the comment-lines from the input Hotspots Table and pass it along
    # to the output Hotspots Table.
    comment_lines = [line for line in open(options.in_table, 'rU').readlines() if line.startswith('#')]

    write_hotspot_table(options.out_table, hotspots_df, comment_lines)

def import_hotspots_table(infile):
    """Import and return the Hotspots Table as a panda dataframe."""

    in_df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    in_comment_lines = [line for line in open(infile, 'rU').readlines() if line.startswith('#')]

    return in_df, in_comment_lines

def write_hotspot_table(outfile, df, comment_lines):
    """Write the Hotspots Table to file."""

    fout = open(outfile, 'w')

    for line in comment_lines:
        fout.write(line)

    df.to_csv(fout, sep="\t", na_rep='nan', index=False)
    fout.close()

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to input Mutation Count file.")

    parser.add_argument("--out_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to output Mutation Count file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_table:", repr(options.in_table)
    print "##   out_table:", repr(options.out_table)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
