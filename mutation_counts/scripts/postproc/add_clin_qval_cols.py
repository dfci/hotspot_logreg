#!/usr/bin/env python2

"""
SYNOPSIS

    Compute and add Clinical QValue columns to the Clinical Hotspot Table.

NOTES

    Since (A) there are multiple hotspots and (B) we are performing multiple
    clinical association (enrichment) tests on each individual hotspot, will
    need to perform multiple hypothesis testing correction on the reported
    p-values.

    We are currently computing the adjusted p-values (q-values) using Benjamini
    and Hochberg (1995) method, which assumes independence or positive
    correlation between the p-values.

    As a starting point, the correction is performed seperately for each clinical
    association test (e.g. ER+ vs. ER-, HER2+ vs. HER2- and etc). So we are
    currently:
        (A) Correcting for the fact that there are multiple hotspots
        (B) But not correcting for the the fact that multiple association tests
            are performed on each individual hotspot.

EXAMPLES

    # (1) Compute and add Clinical QValue columns for Private + Public Analysis
    ./add_clin_qval_cols.py \
        --in_table ../private/reliance_point_outputs/05May2017/05May2017_ReliancePoint_ClinicalHotspotsTable.CompareTaylor.txt \
        --out_table tmp/05May2017_ReliancePoint_ClinicalHotspotsTable.CompareTaylor.WithQValue.txt

    # (2) Compute and add Clinical QValue columns for Public-Only Analysis
    ./add_clin_qval_cols.py \
        --in_table ../private/count_tables/04May2017.count_table.merged.only_public.use_hotspot_stats_model.with_breakdown.CompareTaylor.txt \
        --out_table tmp/04May2017.count_table.merged.only_public.use_hotspot_stats_model.with_breakdown.CompareTaylor.WithQValue.txt

EDGE_CASES


AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import argparse
import numpy as np
import pandas as pd
import statsmodels
from statsmodels.sandbox.stats.multicomp import multipletests


def main(options):

    table_df, comment_lines = import_hotspots_table(options.in_table)

    in_colnames = list(table_df)
    pval_colnames = list()
    for colname in in_colnames:
        if not colname.endswith('PValue') or colname == 'Is_Hotspot_PValue':
            continue
        pval_colnames.append(colname)

    out_colnames = in_colnames
    num_hotspots = len(table_df.index)

    # generate a list of valid pValues and track their associated test and row index
    all_pvalues = []
    track_tests_dict = {
        col: {
            'num_rows': 0,
            'idxs': []
        }
        for col in pval_colnames
    }
    for col in pval_colnames:

        print "##\n## INFO: pval_colname:", repr(col)

        # create column in table
        qval_colname = col.replace("PValue", "QValue")
        col_index = out_colnames.index(col)
        out_colnames.insert(col_index + 1, qval_colname)

        # get pValues
        f1 = (table_df[col] != '-')
        pvals = table_df[f1][col].tolist()
        track_tests_dict[col]['idxs'] = table_df[f1].index

        all_pvalues.extend([float(i) for i in pvals])
        track_tests_dict[col]['num_rows'] = len(pvals)

        print '## INFO: Number of hotspots analyzed: %d' % len(pvals)

    # compute qValues
    qvals = compute_q_values(all_pvalues)

    # dissociate qvals to their respective tests and indices (preserve test order)
    for col in pval_colnames:

        qcol = col.replace('PValue', 'QValue')
        all_qvalues = ['-'] * num_hotspots

        if track_tests_dict[col]['num_rows'] != 0:
            this_tests_valid_qvals = qvals[:track_tests_dict[col]['num_rows']]
            i = 0
            for idx in track_tests_dict[col]['idxs']:
                all_qvalues[idx] = this_tests_valid_qvals[i]
                i += 1

        table_df[qcol] = ['%.3e' % i if isinstance(i, float) else i for i in all_qvalues]

    # Ensure that there is no duplicated columns
    assert len(set(out_colnames)) == len(out_colnames)

    # Reorder the columns to be in the desired order.
    table_df = table_df[out_colnames]

    # Write the output Hotspots Table
    write_hotspot_table(options.out_table, table_df, comment_lines)


def compute_q_values(pvalues):
    """
    Compute the FDR-adjusted p-values (q-values) from the input list of
    uncorrected p-values.

    Return the q-values as a numpy array.

    Note
    ----
    (1) Take care to maintain the order of the elements in the list.
    """

    # Convert the p-values from string to float
    for index in range(len(pvalues)):
        pvalue = pvalues[index]

        assert pvalue != 'inf'
        assert pvalue != 'nan'

        try:
            pvalue = float(pvalue)
        except:
            raise Exception("Unable to convert pvalue (%s) to a float" % pvalue)

        assert pvalue > 0.0
        assert pvalue <= 1.0

        pvalues[index] = pvalue

    pvalues = np.asarray(pvalues)

    # fdr_bh: Benjamini/Hochberg  (non-negative)
    _, qvalues, _, _  = statsmodels.sandbox.stats.multicomp.multipletests(pvalues,
                                                                          method='fdr_bh',
                                                                          is_sorted=False)

    return qvalues


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
                        help="Path to the input hotspots table.")

    parser.add_argument("--out_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output hotpots table.")

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
