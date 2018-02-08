#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the Sanger Institute WGS Breast Cancer Clinical data into a format
    that is suitable for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

EXAMPLES

    ./prep_sanger_wgs_clinical.py \
        --in_clinical ../public/clinical/raw/sanger_wgs/nature17676-s3/Supplementary_Table_1.CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv \
        --in_maf ../public/from_intelccc/sanger_somatic_coding_annotated.maf \
        --out_clinical ../public/clinical/sanger_wgs.clinical_data.tsv

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Feb-2017)
"""


import sys
import os
import time
import copy
import argparse

import numpy as np
import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

OUT_COLS = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status', 'HER2_Status', 'Biopsy_Site_Type', 'Histology_Type', 'Gender']

RECEPTOR_STATUS_COLS = ['ER_Status', 'PR_Status', 'HER2_Status']
RECEPTOR_STATUSES = ['Positive', 'Negative', 'Unknown']

BIOPSY_SITE_TYPES = ['Primary', 'Local_Recurrence', 'Metastatic', 'Unknown']

HISTOLOGY_TYPES = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                   'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

GENDER_TYPES = ['Female', 'Male', 'Unknown']

def main(options):

    clinical_df = import_sanger_wgs_clinical(options.in_clinical)

    maf_sample_barcodes = import_maf_sample_barcodes(options.in_maf)

    suffix_counts = dict()

    clinical_df['Tumor_Sample_Barcode'] = None

    for index, row in clinical_df.iterrows():

        sample_id = row['Sample_ID']

        match_sample_barcodes = list()

        # Find the sample_barcode that has sample_id as its prefix.
        for sample_barcode in maf_sample_barcodes:
            if sample_barcode.startswith(sample_id):
                match_sample_barcodes.append(sample_barcode)

        if len(match_sample_barcodes) == 1:
            match_sample_barcode = match_sample_barcodes[0]
            suffix = match_sample_barcode[len(sample_id):]

            if suffix not in suffix_counts:
                suffix_counts[suffix] = 0
            suffix_counts[suffix] += 1

            # Update the Sample_ID to match the MAF's Tumor_Sample_Barcode.
            clinical_df.set_value(index, 'Tumor_Sample_Barcode', match_sample_barcode)

        elif len(match_sample_barcodes) == 0:
            # Consistency check. Ensure that there is 1-to-1 mapping from Clinical's
            # sample_id and MAF's Tumor_Sample_Barcode.
            raise Exception("Cannot find matching MAF tumor_sample_barcode for clinical sample_id ('%s')." % sample_id)

        else: #  len(match_sample_barcodes) > 1:
            # Consistency check. Ensure that there is 1-to-1 mapping from Clinical's
            # sample_id and MAF's Tumor_Sample_Barcode.
            raise Exception("Mutiple matching MAF tumor_sample_barcodes for clinical sample_id ('%s')." % sample_id)

    # Add Center information
    clinical_df['Center'] = 'SangerWGS'

    # Keep only the desired output columns in the desired order.
    clinical_df = clinical_df[OUT_COLS]

    # Ensure that there is no missing data in any of the columns.
    assert not clinical_df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not clinical_df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Ensure that every sample_barcode in the MAF file maps to a clinical sample_id
    assert set(clinical_df['Tumor_Sample_Barcode']) == set(maf_sample_barcodes)

    # Output the Clinical Data
    clinical_df.to_csv(options.out_clinical, sep="\t", na_rep='', index=False)

    print "##", "-" * 50
    print "## Outfile Summary"
    print "##   Num Samples Clinical Total: %s" % len(clinical_df)
    print "##   Num Samples Genomics Total: %s" % len(maf_sample_barcodes)
    print "##   Suffix Counts Breakdown: %s" % suffix_counts
    print "##", "-" * 50

def import_sanger_wgs_clinical(infile):
    """Import the Sanger Institute WGS clinical data from Supplemental Table 1
    of the paper.

    Notes
    -----
    - Extract the following data in standardized format:
         (1)  Sample_ID (e.g. 'PD10010')
         (2)  ER_Status (see RECEPTOR_STATUSES list)
         (3)  PR_Status (see RECEPTOR_STATUSES list)
         (4)  HER2_Status (see RECEPTOR_STATUSES list)
         (5)  Biopsy_Site_Type (see BIOPSY_SITE_TYPES list)
         (6)  Histology_Type (see HISTOLOGY_TYPES list)
         (7)  Gender (see GENDER_TYPES list)
    """

    df = pd.read_csv(infile, skiprows=1)

    # Check that all the require columns exist
    required_columns = ['sample_name', 'final.ER', 'final.PR', 'final.HER2', 'specimen_type', 'Histopathological_subtype', 'donor_gender']

    df = df[required_columns]

    # Standardize the column names
    df.rename(columns={'sample_name': 'Sample_ID'}, inplace=True)
    df.rename(columns={'final.ER': 'ER_Status'}, inplace=True)
    df.rename(columns={'final.PR': 'PR_Status'}, inplace=True)
    df.rename(columns={'final.HER2': 'HER2_Status'}, inplace=True)
    df.rename(columns={'specimen_type': 'Biopsy_Site_Type'}, inplace=True)
    df.rename(columns={'Histopathological_subtype': 'Histology_Type'}, inplace=True)
    df.rename(columns={'donor_gender': 'Gender'}, inplace=True)

    # Fill missing ER/PR/HER2 status values
    for colname in RECEPTOR_STATUS_COLS:
        df[[colname]] = df[[colname]].fillna(value='Unknown')

    # Standardize the values in the [ER,PR,HER2]_Receptor_Status columns
    for colname in RECEPTOR_STATUS_COLS:
        df[colname].replace('positive', 'Positive', inplace=True)
        df[colname].replace('negative', 'Negative', inplace=True)
        df[colname].replace('no_data_supplied', 'Unknown', inplace=True)

        # Ensure that there is no unexpected receptor_status values.
        assert set(df[colname].unique()) <= set(RECEPTOR_STATUSES)

    # Standardize the values in the 'Biopsy_Site_type' column
    df['Biopsy_Site_Type'].replace('tumour_primary', 'Primary', inplace=True)
    df['Biopsy_Site_Type'].replace('tumour_local_recurrence', 'Local_Recurrence', inplace=True)
    df['Biopsy_Site_Type'].replace('metastasis', 'Metastatic', inplace=True)

    # Standardize the values in the 'Histology_Type' column
    #     (1) Used the following 3 sources to guide the mapping:
    #           (A) http://emedicine.medscape.com/article/1954658-overview
    #           (B) http://breast-cancer.ca/  (navigate to Section 5)
    #           (C) http://www.breastcancer.org/symptoms/types
    #
    #     (2) Assume that all tumors are NOT In Situ Carcinoma.
    df['Histology_Type'].replace('ductal',                   'Invasive_Ductal_Carcinoma',        inplace=True)
    df['Histology_Type'].replace('lobular',                  'Invasive_Lobular_Carcinoma',       inplace=True)
    df['Histology_Type'].replace('adenoid_cystic_carcinoma', 'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('apocrine',                 'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('cribriform_tubular',       'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('lobular_Pleomorphic',      'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('micropapillary',           'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('mucinous',                 'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('neuroendocrine',           'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('papillary',                'Other_Invasive_Breast_Carcinoma',  inplace=True)
    df['Histology_Type'].replace('metaplastic',              'Other_Breast_Cancer',              inplace=True)
    df['Histology_Type'].replace('no_data_supplied',         'Unknown',                          inplace=True)

    # Standardize the values in the 'Gender column
    df['Gender'].replace('female', 'Female', inplace=True)
    df['Gender'].replace('male', 'Male', inplace=True)

    # Ensure that there is no unexpected 'Biopsy_Site_Type' values.
    assert set(df['Biopsy_Site_Type'].unique()) <= set(BIOPSY_SITE_TYPES)

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df['Histology_Type'].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there is no unexpected 'Gender' values.
    assert set(df['Gender'].unique()) <= set(GENDER_TYPES)

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Sample_ID' values
    assert not df.duplicated(subset=['Sample_ID']).any()

    return df

def import_maf_sample_barcodes(infile):
    """Import a set of tumor sample barcodes from the MAF file.
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header = 0)

    sample_barcodes = list(df['Tumor_Sample_Barcode'].unique())

    return sample_barcodes

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Clinical Data.")

    parser.add_argument("--in_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to input MAF file.")

    parser.add_argument("--out_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to output Clinical Data.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_clinical: ", repr(options.in_clinical)
    print "##   in_maf: ", repr(options.in_maf)
    print "##   out_clinical: ", repr(options.out_clinical)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
