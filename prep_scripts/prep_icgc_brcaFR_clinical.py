#!/usr/bin/env python

"""
SYNOPSIS
    Prepare the ICGC BRCA-FR donor and sample files into a clinical tsv format that is suitable for input into
    the Hotspots 'Mutation_Counts' pipeline.

NOTES
    - From https://dcc.icgc.org/projects/BRCA-FR:
        Tumour Subtype:	Subtype defined by an amplification of the HER2 gene
        so we can assing HER2+ to this study's clinical file.

    - Data was downloaded from https://dcc.icgc.org/releases/current/Projects/BRCA-FR.

EXAMPLES
    ```
    ./prep_icgc_brcaFR_clinical.py \
        --donor-file ./donor.BRCA-FR.tsv \
        --sample-file ./sample.BRCA-FR.tsv \
        --out-clinical ./icgc_brca_fr.clinical.r1.tsv
    ```

AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import json
import argparse

import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

clinical_columns = [
    'Tumor_Sample_Barcode',
    'Center',
    'ER_Status',
    'PR_Status',
    'HER2_Status',
    'Biopsy_Site_Type',
    'Histology_Type',
    'Gender'
]

donor_columns = [
    'icgc_donor_id',
    'project_code',
    'study_donor_involved_in ',
    'submitted_donor_id',
    'donor_sex',
    'donor_vital_status',
    'disease_status_last_followup',
    'donor_relapse_type',
    'donor_age_at_diagnosis',
    'donor_age_at_enrollment',
    'donor_age_at_last_followup',
    'donor_relapse_interval',
    'donor_diagnosis_icd10',
    'donor_tumour_staging_system_at_diagnosis',
    'donor_tumour_stage_at_diagnosis',
    'donor_tumour_stage_at_diagnosis_supplemental',
    'donor_survival_time',
    'donor_interval_of_last_followup',
    'prior_malignancy',
    'cancer_type_prior_malignancy',
    'cancer_history_first_degree_relative'
]

sample_columns = [
    'icgc_sample_id',
    'project_code',
    'submitted_sample_id',
    'icgc_specimen_id',
    'submitted_specimen_id',
    'icgc_donor_id',
    'submitted_donor_id',
    'analyzed_sample_interval',
    'percentage_cellularity',
    'level_of_cellularity',
    'study',
]

biopsy_site_dict = {
    'distant recurrence/metastasis': 'Metastatic',
    'local recurrence': 'Primary'
}


class Clinical:
    """
    For this project, we have no hormone receptor status or histology
    type information.

    Gender was pulled from the donor file

    Required columns are:
     - Tumor_Sample_Barcode
     - Center
     - ER_Status
     - PR_Status
     - HER2_Status
     - Biopsy_Site_Type
     - Histology_Type
     - Gender
    """

    def __init__(self, args):
        self.args = args

        self.donor_df = None
        self.sample_df = None
        self.clinical_df = None

        # iterator used for logging
        self.i = 0

    def get_gender(self, donor_id):
        """
        Returns gender from icgc_sample_id

        :param donor_id: ICGC Sample ID
        :return: Male or Female
        """
        f1 = (self.donor_df.icgc_donor_id == donor_id)
        return self.donor_df[f1]['donor_sex'].tolist()[0].title()

    @staticmethod
    def get_biopsy_site(donor_relapse_type):
        """
        Convert the values in the 'donor_relapse_type' to biopsy site type values that
        will be accepted by the Hotspot 'Mutation_Counts' pipeline.
        """

        if pd.isnull(donor_relapse_type):
            return 'Unknown'
        else:
            return biopsy_site_dict[donor_relapse_type]

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'

        # Donor file
        if self.args.donor_file:
            if not self.args.sample_file:
                print '## ERROR: To generate the clinical file, you must provide the donor and sample '
                print '##        file paths because the donor ID is mapped to the sample ID in the'
                print '##        mutations file through these.'

            self.clinical_df = pd.DataFrame(columns=clinical_columns)
            self.donor_df = pd.read_csv(self.args.donor_file, sep='\t', header=0)
            self.sample_df = pd.read_csv(self.args.sample_file, sep='\t', header=0)

    def write_results(self):
        """Write clinical dataframe to file."""

        # clinical file
        with open(self.args.out_clinical, 'w') as ff:
            ff.write('#version 1.0\n')
            self.clinical_df.to_csv(ff, sep='\t', index=False)

    def create_clinical_tsv(self):
        """
        Creates clinical file from the three dataframes donor_df, specimen_df, and sample_df,
        as well as assumes hormone receptor status and histology type through the title.
        """

        print '\n## Creating clinical file...'
        print'## Sample file contains %d samples' % len(self.sample_df.index)
        self.clinical_df.Tumor_Sample_Barcode = self.sample_df.icgc_sample_id
        self.clinical_df.Center = self.sample_df.project_code
        self.clinical_df.ER_Status = 'Unknown'
        self.clinical_df.PR_Status = 'Unknown'
        self.clinical_df.HER2_Status = 'Positive'
        self.clinical_df.Biopsy_Site_Type = self.donor_df.donor_relapse_type.apply(lambda x: self.get_biopsy_site(x))
        self.clinical_df.Histology_Type = 'Unknown'
        self.clinical_df.Gender = self.sample_df.icgc_donor_id.apply(lambda x: self.get_gender(x))
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

    parser.add_argument('--donor-file', dest='donor_file', required=False,
                        help='Path to donor clinical file.')

    parser.add_argument('--sample-file', dest='sample_file', required=False,
                        help='Path to sample clinical file.')

    parser.add_argument('--out-clinical', dest='out_clinical', required=False,
                        help='Path to output clinical file.')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
