#!/usr/bin/env python2

"""
SYNOPSIS
    This script splits the GENIE .bed file found here:
        https://www.synapse.org/#!Synapse:syn7844527

EXAMPLES

    ./prep_genie_bed.py \
        --in-bed ./genie_combined.bed \
        --out-dir .

AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (June-2017)
"""

import sys
import os
import json
import time
import argparse
import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)


class ParseBed:
    def __init__(self, args):
        self.args = args
        self.genie_bed_df = None
        self.genie_dict = {
            'vicc': pd.DataFrame(),
            'grcc': pd.DataFrame(),
            'uhn': pd.DataFrame(),
            'mda': pd.DataFrame()
        }

    def load_bed(self):
        self.genie_bed_df = pd.read_csv(self.args.in_bed, sep='\t')

    def split_bed(self):

        for center in self.genie_dict:
            f1 = (self.genie_bed_df.SEQ_ASSAY_ID.str.startswith(center.upper()))
            self.genie_dict[center] = self.genie_bed_df[f1]

    def write_results(self):

        for center in self.genie_dict:
            filename = '%s/GENIE.%s.r1.bed' % (self.args.out_dir, center)
            with open(filename, 'w') as ff:
                ff.write('#version 1.0\n#')
                self.genie_dict[center].to_csv(ff, sep='\t', index=False)


def main(args):
    run = ParseBed(args)
    run.load_bed()
    run.split_bed()
    run.write_results()


if __name__ == '__main__':
    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-bed', dest='in_bed', required=True,
                        help='Path to GENIE .bed file')

    parser.add_argument('--out-dir', dest='out_dir', required=True,
                        help='Path to output directory for GRCC, VICC, MDA and UHN .bed files')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
