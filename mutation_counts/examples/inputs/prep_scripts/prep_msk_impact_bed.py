#!/usr/bin/env python2

"""
SYNOPSIS
    Extract MSK-IMPACT panels' target interval content from the combined GENIE
    BED file and preprocess the content into format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (a) Extract the bed regions by matching panel name value in the 'SEQ_ASSAY_ID'
        column.

    (b) Prepend '#' character to the header line.

EXAMPLES

    # Create MSK-IMPACT341 Target BED file
    ./prep_msk_impact_bed.py \
        --in_bed ../public/genie/raw/03Feb2017/genie_combined.bed \
        --seq_assay_id MSK-IMPACT341 \
        --out_bed ../public/genie/genie_msk.IMPACT341.r1.check.bed

    # Create MSK-IMPACT410 Target BED file
    ./prep_msk_impact_bed.py \
        --in_bed ../public/genie/raw/03Feb2017/genie_combined.bed \
        --seq_assay_id MSK-IMPACT410 \
        --out_bed ../public/genie/genie_msk.IMPACT410.r1.check.bed


AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (April-2017)
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
            # Input Column Header line.
            assert cols == ["Chromosome", "Start_Position", "End_Position", "Hugo_Symbol",
                            "SEQ_ASSAY_ID", "Feature_Type"]

            cols[0] = "#" + cols[0]

        else:
            # Data line.
            chrom = cols[0]
            seq_assay_id = cols[4]

            assert chrom in set(['%s' % chrom for chrom in range(1,23)] + ['X', 'Y', 'MT'])

            if seq_assay_id != options.seq_assay_id:
                continue

        fout.write("\t".join(cols) + "\n")

    fout.close()

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to input BED file")

    parser.add_argument("--seq_assay_id", action="store", required=True,
                         help="Filter for only data of specified Panel.")

    parser.add_argument("--out_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output BED file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_bed: ", repr(options.in_bed)
    print "##   seq_assay_id: ", repr(options.seq_assay_id)
    print "##   out_bed: ", repr(options.out_bed)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
