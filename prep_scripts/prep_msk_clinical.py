#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the MSK 10k data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

    As of June 12th, 2017, the cBioPortal team did not include the SEQ_ASSAY_ID column in the data
    uploaded to the cBioPortal Datahub, so in order to map the panel version in the 10K dataset I parsed
    the last 3 digits of the SAMPLE_ID as follows:
        - 'IM3' == 'IMPACT341'
        - 'IM5' == 'IMPACT410'

NOTES

    - This study has no hormone receptor status

    - Public dataset downloaded from cBioPortal datahub:
    https://github.com/cBioPortal/datahub/blob/master/public/msk_impact_2017.tar.gz

EXAMPLES

    ./prep_msk_clinical.py \
        --in-sample ./data_clinical_sample.txt \
        --in-patient ./data_clinical_patient.txt \
        --out-dir .

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

biopsy_site_type_dict = {
    'Primary': 'Primary',
    'Metastasis': 'Metastatic'
}

histology_type_dict = {
    'IDC': 'Invasive_Ductal_Carcinoma',
    'ILC': 'Invasive_Lobular_Carcinoma',
    'MDLC': 'Mixed_Ductal_and_Lobular_Carcinoma',
    'Breast Carcinoma': 'Other_Invasive_Breast_Carcinoma',
    'Breast Sarcoma': 'Other_Breast_Cancer'
}


class Clinical:

    def __init__(self, args):

        self.args = args
        self.sample_df = None
        self.patient_df = None
        self.impact341_df = pd.DataFrame(columns=OUT_COLS)
        self.impact410_df = pd.DataFrame(columns=OUT_COLS)

        self.panel_dict = {
            'MSK-IMPACT341': self.impact341_df,
            'MSK-IMPACT410': self.impact410_df
        }

        self.oncotree = oncotreenx.build_oncotree()
        self.i = 0

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.sample_df = pd.read_csv(self.args.in_sample, sep='\t', skiprows=4)
        self.patient_df = pd.read_csv(self.args.in_patient, sep='\t', skiprows=4)

    def write_results(self, df, center):
        """Write clinical dataframe to file."""

        # clinical file
        filename = '%s/msk10k.%s.clinical.r2.tsv' % (self.args.out_dir, center.replace('MSK-', '').lower())
        with open(filename, 'w') as ff:
            ff.write('#version 1.0\n')
            df.to_csv(ff, sep='\t', index=False)

    def get_gender(self, sample_id):
        """Return gender"""

        # for logging
        self.i += 1
        if self.i % 2500 == 0:
            print '## %s rows processed' % self.i

        f1 = (self.sample_df.SAMPLE_ID == sample_id)
        patient_id = self.sample_df[f1]['PATIENT_ID'].tolist()[0]

        f2 = (self.patient_df.PATIENT_ID == patient_id)
        return self.patient_df[f2]['SEX'].tolist()[0]

    @staticmethod
    def get_biopsy_site_type(sample_type):
        """Return the corresponding biopsy type value that the hotspot pipeline accepts"""

        if pd.isnull(sample_type):
            return 'Unknown'
        else:
            return biopsy_site_type_dict[sample_type]

    def get_histology_type(self, oncotree_code):
        """Convert oncotree code to histology type"""

        if oncotree_code in histology_type_dict:
            return histology_type_dict[oncotree_code]

        try:
            metamaintype = self.oncotree.node[oncotree_code]['metamaintype']
        except KeyError:
            print '## WARNING: \'%s\' is not a valid oncotree code. Sample was removed from analysis' % oncotree_code
            return 'removeme'

        if metamaintype in histology_type_dict:
            return histology_type_dict[metamaintype]
        else:
            return 'removeme'

    def split_by_panel(self):

        f_341 = (self.sample_df.SAMPLE_ID.str.endswith('IM3'))
        f_410 = (self.sample_df.SAMPLE_ID.str.endswith('IM5'))

        self.panel_dict['MSK-IMPACT341'] = self.sample_df[f_341]
        self.panel_dict['MSK-IMPACT410'] = self.sample_df[f_410]

    def create_clinical_tsv(self):

        print '\n## Creating clinical file...'

        for center, df in self.panel_dict.iteritems():
            print'## Sample file contains %d samples' % len(df.SAMPLE_ID.unique())
            clinical_df = pd.DataFrame(columns=OUT_COLS)
            clinical_df.Tumor_Sample_Barcode = df.SAMPLE_ID
            clinical_df.Center = center
            clinical_df.ER_Status = 'Unknown'
            clinical_df.PR_Status = 'Unknown'
            clinical_df.HER2_Status = 'Unknown'
            clinical_df.Biopsy_Site_Type = df.SAMPLE_TYPE.apply(lambda x: self.get_biopsy_site_type(x))
            clinical_df.Histology_Type = df.ONCOTREE_CODE.apply(lambda x: self.get_histology_type(x))
            clinical_df.Gender = df.SAMPLE_ID.apply(lambda x: self.get_gender(x))
            clinical_df = clinical_df.fillna(value='Unknown', axis=1)

            # drop non breast cancers
            f1 = (clinical_df.Histology_Type != 'removeme')
            clinical_df = clinical_df[f1]

            self.write_results(clinical_df, center)


def main(args):

    run = Clinical(args)
    run.load_data()
    run.split_by_panel()
    run.create_clinical_tsv()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-sample', dest='in_sample', required=True,
                        help='Path to input MSK sample file')

    parser.add_argument('--in-patient', dest='in_patient', required=True,
                        help='Path to input MSK patient file')

    parser.add_argument('--out-dir', dest='out_dir', required=True,
                        help='Path to output clinical file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
