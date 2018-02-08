#!/usr/bin/env python2

"""
SYNOPSIS
     Prepares the GENIE data into a format that is suitable for input
     into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    - This study has no hormone receptor status

    - Public dataset downloaded from Synapse:
    https://www.synapse.org/#!Synapse:syn7851246

EXAMPLES
    ```
    ./prep_genie_clinical.py \
        --in-clinical ./data_clinical.txt \
        --out-clinical .
    ```


AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import json
import time
import argparse
import pandas as pd

import oncotreenx

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

allowed_cancer_types = ['Breast Carcinoma', 'Breast Sarcoma']

histology_type_dict = {
    'Breast Invasive Lobular Carcinoma': 'Invasive_Lobular_Carcinoma',
    'Breast Invasive Ductal Carcinoma': 'Invasive_Ductal_Carcinoma',
    'Breast Mixed Ductal and Lobular Carcinoma': 'Mixed_Ductal_and_Lobular_Carcinoma',
    'Breast Carcinoma': 'Other_Invasive_Breast_Carcinoma',
    'Breast Sarcoma': 'Other_Breast_Cancer'
}

biopsy_site_type_dict = {
    'Primary': 'Primary',
    'Metastasis': 'Metastatic'
}


class Clinical:

    def __init__(self, args):

        self.args = args

        self.mda_df = None
        self.vicc_df = None
        self.uhn_df = None
        self.grcc_df = None

        self.oncotree = oncotreenx.build_oncotree()

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.clinical_df = pd.read_csv(self.args.in_clinical, sep='\t')

    @staticmethod
    def write_results(df, center):
        """Write clinical dataframe to file."""

        # clinical file
        filename = '%s/%s.clinical.r1.tsv' % (args.out_clinical, center)
        with open(filename, 'w') as ff:
            ff.write('#version 1.0\n')
            df.to_csv(ff, sep='\t', index=False)

    def get_histology_type(self, cancer_type_text):
        """Convert oncotree code to histology type"""

        if cancer_type_text in histology_type_dict:
            return histology_type_dict[cancer_type_text]

        try:
            metamaintype = self.oncotree.node[oncotreenx.lookup_text(self.oncotree, cancer_type_text)]['metamaintype']
        except KeyError:
            print '## WARNING: \'%s\' is not a valid oncotree text input.' \
                  ' Sample was removed from analysis' % cancer_type_text
            return 'removeme'

        if metamaintype in histology_type_dict:
            return histology_type_dict[metamaintype]
        else:
            return 'removeme'

    def split_by_center(self):
        """
        Split clinical dataframe into the following centers:

        - MDA   | University of Texas, MD Anderson Cancer Center
        - VICC  | Vanderbilt-Ingram Cancer Center
        - UHN   | Princes Margaret Cancer Center, University Health Network
        - GRCC  | Institute Gustave Roussy, France
        """

        f1 = (self.clinical_df.CENTER == 'MDA')
        self.mda_df = self.clinical_df[f1]
        print '## INFO: MDA has %s samples' % len(self.mda_df.index)

        f2 = (self.clinical_df.CENTER == 'VICC')
        self.vicc_df = self.clinical_df[f2]
        print '## INFO: VICC has %s samples' % len(self.vicc_df.index)

        f3 = (self.clinical_df.CENTER == 'UHN')
        self.uhn_df = self.clinical_df[f3]
        print '## INFO: UHN has %s samples' % len(self.uhn_df.index)

        f4 = (self.clinical_df.CENTER == 'GRCC')
        self.grcc_df = self.clinical_df[f4]
        print '## INFO: GRCC has %s samples' % len(self.grcc_df.index)

    @staticmethod
    def get_biopsy_site_type(sample_type):
        """Return the corresponding biopsy type value that the hotspot pipeline accepts"""

        if pd.isnull(sample_type) or sample_type == 'Unspecified':
            return 'Unknown'
        else:
            return biopsy_site_type_dict[sample_type]

    def create_clinical_tsv(self):
        """Create the clinical tsv for each center and write the results"""

        for df in [self.mda_df, self.vicc_df, self.uhn_df, self.grcc_df]:

            clinical_df = pd.DataFrame(columns=OUT_COLS)

            center = df.CENTER.unique().tolist()[0]
            print '\n## Creating clinical file %s...' % center
            print'## Sample file contains %d samples' % len(df.SAMPLE_ID.unique())
            clinical_df.Tumor_Sample_Barcode = df.SAMPLE_ID
            clinical_df.Center = df.CENTER
            clinical_df.ER_Status = 'Unknown'
            clinical_df.PR_Status = 'Unknown'
            clinical_df.HER2_Status = 'Unknown'
            clinical_df.Biopsy_Site_Type = df.SAMPLE_TYPE.apply(lambda x: self.get_biopsy_site_type(x))
            clinical_df.Histology_Type = df.CANCER_TYPE_DETAILED.apply(lambda x: self.get_histology_type(x))
            clinical_df.Gender = 'Unknown'
            clinical_df = clinical_df.fillna(value='Unknown', axis=1)

            # drop non breast cancers
            clinical_df = clinical_df[clinical_df.Histology_Type != 'removeme']

            self.write_results(clinical_df, center)


def main(args):

    run = Clinical(args)
    run.load_data()
    run.split_by_center()
    run.create_clinical_tsv()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-clinical', dest='in_clinical', required=True,
                        help='Path to input GENIE clinical file')

    parser.add_argument('--out-clinical', dest='out_clinical', required=True,
                        help='Path to output GENIE clinical files. Please input'
                             'a path and not a filename, as this script will generate'
                             'four files for each center.')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
