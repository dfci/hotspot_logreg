#!/usr/bin/env python2

"""
SYNOPSIS

    Annotate each putative hotspot in the Private+Public Table with information
    from the Public-Only Table.

NOTES

    (1) Add the following annotation 3 columns:

        (A) "PublicOnly_Num_Samples_with_Mutation":
                Number of Samples with Mutation from Public Data Only Analysis Table.

        (B) "PublicOnly_Is_Hotspot_PValue":
                Pvalue of the Hotspot from the Public Data Only Analysis Table.

        (C) "PublicOnly_Is_Hotspot_QValue":
                QValue of the Hotspot from the Public Data Only Analysis Table.

    (2) Match the putative hotspots in the two tables by (Gene_Symbol, Transcript_ID, Ref_Codon)
        tuple key.

EXAMPLES

    # (1)
    ./add_public_only_cols.py \
        --public_only_table ../private/count_tables/05May2017.count_table.merged.only_public.with_breakdown.no_qvalue_filter.min_count_1.txt \
        --in_table ../private/reliance_point_outputs/05May2017/05May2017_ReliancePoint_ClinicalHotspotsTable.CompareTaylor.WithQValue.txt \
        --out_table tmp/05May2017_ReliancePoint_ClinicalHotspotsTable.CompareTaylor.WithQValue.ComparePublicOnly.txt

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import argparse
import pandas as pd

def main(options):

    public_only_df, _ = import_hotspots_table(options.public_only_table)
    full_df, comment_lines = import_hotspots_table(options.in_table)

    public_only_df.index = (public_only_df['Gene_Symbol'] +
                            "_" + public_only_df['Transcript_ID'] +
                            "_" + public_only_df['Ref_Codon'])

    colnames = ['Num_Samples_with_Mutation', 'Is_Hotspot_PValue', 'Is_Hotspot_QValue']
    public_only_df = public_only_df[colnames]

    for colname in colnames:
        public_only_df.rename(columns={colname: 'PublicOnly_' + colname}, inplace=True)

    full_df.index = full_df['Gene_Symbol'] + "_" + full_df['Transcript_ID'] + "_" + full_df['Ref_Codon']

    # Perform left outer join.
    full_df = pd.merge(full_df, public_only_df,
                       how="left",  # LEFT OUTER JOIN
                       left_index=True,
                       right_index=True,
                       sort=False)


    # Write the Output Hotspots Table
    write_hotspot_table(options.out_table, full_df, comment_lines)

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

    parser.add_argument("--public_only_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to the public_only_table")

    parser.add_argument("--in_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input mutation count table (transfer data to)")

    parser.add_argument("--out_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output mutation count table.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   public_only_table:", repr(options.public_only_table)
    print "##   in_table:", repr(options.in_table)
    print "##   out_table:", repr(options.out_table)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
