#!/usr/bin/env python2

"""
SYNOPSIS
    Preprocess the DFCI-ONCOPANEL-2 BED file into a format that is suitable for
    input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (a) Remove the 'chr' prefix from the 'chrom' (1st) column

    (b) Add '#' character to the column header line.

EXAMPLES

    ./prep_dfci_oncopanel2_bed.py \
        --in_bed ../private/raw/29Nov2016/PROFILE_POPv2.filtered.bed \
        --out_bed ../private/dfci.PROFILE_POPv2.filtered.bed

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Dec-2016)
"""


import sys
import os
import time
import copy
import subprocess
import argparse

import pandas as pd

def main(options):

    fout = open(options.out_bed, 'w')

    line_num = 0
    for line in open(options.in_bed):
        line_num += 1

        line = line.rstrip('\r\n')
        cols = line.split('\t')

        if line_num == 1:
            # Column Header line.
            assert cols == ["chrom","chromStart","chromStop","name"]

            cols[0] = "#chrom"

        else:
            # Data line.
            assert cols[0].startswith("chr")
            cols[0] = cols[0][len("chr"):]
            assert cols[0] in set(['%s' % chrom for chrom in range(1,23)] + ['X', 'Y', 'MT'])

        fout.write("\t".join(cols) + "\n")

    fout.close()

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to input BED file")

    parser.add_argument("--out_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output BED file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_bed: ", repr(options.in_bed)
    print "##   out_bed: ", repr(options.out_bed)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
