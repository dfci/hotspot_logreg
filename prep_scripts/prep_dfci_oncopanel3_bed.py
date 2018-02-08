#!/usr/bin/env python2

"""
SYNOPSIS
    Preprocess the DFCI-ONCOPANEL-3 BED file into a format that is suitable for
    input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (a) Skip at meta-data lines

    (b) Assume that start_coord (2nd column) is already in BED format (i.e. 0-based)

    (c) Assume that end_coord (3rd column) is already in BED format (i.e. 1-based)

EXAMPLES

    ./prep_dfci_oncopanel3_bed.py \
        --in_interval_file ../private/raw/06Feb2017/POPv3_AllTargetedRegions.interval_list \
        --out_bed ../private/dfci.PROFILE_POPv3.r1.bed

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Feb-2017)
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

    for line in open(options.in_interval_file):

        line = line.rstrip('\r\n')
        cols = line.split('\t')

        if line.startswith("@"):
            # Skip meta-data line, e.g.
            #       @SQ     SN:22   LN:51304566 .......
            continue

        # Data line
        assert len(cols) > 0
        chrom = cols[0]

        assert chrom in set(['%s' % chrom for chrom in range(1,23)] + ['X', 'Y', 'MT'])


        fout.write("\t".join(cols) + "\n")

    fout.close()

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_interval_file", action="store", required=True,
                        metavar='FILE',
                        help="Path to input Interval file")

    parser.add_argument("--out_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output BED file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_interval_file:", repr(options.in_interval_file)
    print "##   out_bed:", repr(options.out_bed)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
