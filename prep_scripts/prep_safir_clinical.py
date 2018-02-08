#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the SAFIR data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

   - Study is a whole-exome mutational profiling of metastatic breast cancer.
        - 216 samples
        - All metastatic
        - HR/HER2 status for all samples
        - Histology Type unknown for all samples
        - Gender unknown for all samples

   - Paper: https://www.ncbi.nlm.nih.gov/pubmed/28027327
   - Public dataset downloaded from cBioPortal datahub.

EXAMPLES

    ./prep_safir_clinical.py \
        --in-sample ./data_clinical_sample.txt \
        --out-clinical ./safir.clinical_data.tsv

AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (May-2017)
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

OUT_COLS = [
    'Tumor_Sample_Barcode',
    'Center',
    'ER_Status',
    'PR_Status',
    'HER2_Status',
    'Biopsy_Site_Type',
    'Histology_Type',
    'Gender'
]


class Clinical:

    def __init__(self, args):

        self.args = args
        self.sample_df = None
        self.clinical_df = pd.DataFrame(columns=OUT_COLS)

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.sample_df = pd.read_csv(self.args.in_sample, sep='\t', skiprows=4)

    def write_results(self):
        """Write clinical dataframe to file."""

        # clinical file
        with open(self.args.out_clinical, 'w') as ff:
            ff.write('#version 1.0\n')
            self.clinical_df.to_csv(ff, sep='\t', index=False)

    @staticmethod
    def get_her2_status(receptor_status):
        """Return the HER2 status from the IHC_HER2 column in the sample dataframe"""

        if 'her2' not in receptor_status.lower():
            return 'Unknown'
        elif 'her2-' in receptor_status.lower():
            return 'Negative'
        elif 'her2+' in receptor_status.lower():
            return 'Positive'
        else:
            return 'Unknown'

    def create_clinical_tsv(self):
        """
        Creates clinical file from the three dataframes donor_df, specimen_df, and sample_df,
        as well as assumes hormone receptor status and histology type through the title.
        """

        print '\n## Creating clinical file...'
        print'## Sample file contains %d samples' % len(self.sample_df.SAMPLE_ID.unique())
        self.clinical_df.Tumor_Sample_Barcode = self.sample_df.SAMPLE_ID
        self.clinical_df.Center = 'SAFIR'
        self.clinical_df.ER_Status = 'Unknown'
        self.clinical_df.PR_Status = 'Unknown'
        self.clinical_df.HER2_Status = self.sample_df.IHC_HER2.apply(lambda x: self.get_her2_status(x))
        self.clinical_df.Biopsy_Site_Type = 'Metastatic'
        self.clinical_df.Histology_Type = 'Unknown'
        self.clinical_df.Gender = 'Unknown'
        self.clinical_df = self.clinical_df.fillna(value='Unknown', axis=1)


def main(args):

    run = Clinical(args)
    run.load_data()
    run.create_clinical_tsv()
    run.write_results()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-sample', dest='in_sample', required=True,
                        help='Path to input SAFIR sample file')

    parser.add_argument('--out-clinical', dest='out_clinical', required=True,
                        help='Path to output clinical file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
