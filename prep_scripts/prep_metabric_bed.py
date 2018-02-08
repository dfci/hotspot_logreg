#!/usr/bin/env python2

"""
SYNOPSIS
    Preprocess the METABRIC Target BED file into a format that is suitable
    for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Prepend '#' character to header line, so that downstream scripts can correctly
        interpret/ignore this line.

    (2) Remove the 'chr' prefix from the chromosome (1st) column
        (i.e. simple conversion of hg19 to GRCh37).

    (3) Assume all genomics intervals are on standard chromosomes (1,2,3,...,21,22,X,Y,MT)

EXAMPLES

    ./prep_metabric_bed.py \
        --in_bed ../public/raw/metabric.targetedIntervals.txt \
        --out_bed ../public/metabric.targetedIntervals.r1.bed

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Dec-2016)
"""

import sys
import os
import time
import argparse

import pandas as pd

def main(options):

    fout = open(options.out_bed, 'w')

    GRCh37_std_chroms = set(['%s' % chrom for chrom in range(1,23)] + ['X', 'Y', 'MT'])

    line_num = 0
    for line in open(options.in_bed, 'rU'):
        line_num += 1

        # Only '\n' since opened file with universal newline support.
        line = line.rstrip('\n')
        cols = line.split('\t')

        if line_num == 1:
            # Expect column header line.
            assert cols == ["chr", "start", "end", "id", "meanCov"]
            cols[0] = "#chr"

        else:
            # Expect data line.

            # Each data line need to have at least 3 columns ("chr", "start", "end")
            assert len(cols) >= 3

            # Format chromosome column
            assert cols[0].startswith("chr")
            cols[0] = cols[0][len("chr"):]
            assert cols[0] in GRCh37_std_chroms

            # Check that Start Coordinate is a positive integer
            assert cols[1].isdigit()
            assert int(cols[1]) > 0

            # Check that End Coordinate is a positive integer
            assert cols[2].isdigit()
            assert int(cols[2]) > 0

        fout.write('\t'.join(cols) + '\n')

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
