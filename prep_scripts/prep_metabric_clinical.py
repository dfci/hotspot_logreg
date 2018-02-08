#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the METABRIC Breast Cancer Clinical data into a format that is
    suitable for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Support 'details' mode, which will:

          (A) Retain samples which will be normally filtered out (e.g. duplicated
              samples from a patient, Breast Sarcoma samples) and instead
              add 'Exclude_Sample' column  ('True', 'False') and 'Exclude_Reason'
              column (text describing the reason to exclude sample).

              For sample exclusion logic, see filter_clinical_df() function.

        Note
        ----
        Purpose of the 'details' mode is to include extra rows/columns in the
        preprocessed clinical outfile that can be consumed by the 'prep_metabric_maf.py'
        script BUT is not used/supported by the actual mutation counts pipeline.

EXAMPLES

    # (1) Generate the 'regular' preprocessed clinical outfile.
    ./prep_metabric_clinical.py \
        --in_clinical ../public/clinical/raw/metabric/brca_metabric_clinical_data.tsv \
        --out_clinical ../public/clinical/metabric.clinical_data.tsv

    # (2) Generate the 'details' preprocessed clinical outfile.
    ./prep_metabric_clinical.py \
        --in_clinical ../public/clinical/raw/metabric/brca_metabric_clinical_data.tsv \
        --details_mode \
        --out_clinical ../public/clinical/metabric.clinical_data.details.tsv

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Apr-2017)
"""


import sys
import os
import time
import copy
import argparse

import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

OUT_COLS = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status', 'HER2_Status', 'Biopsy_Site_Type', 'Histology_Type']
DETAILS_OUT_COLS = ['Exclude_Sample', 'Exclude_Reason']

BIOPSY_SITE_TYPES = ['Primary', 'Local_Recurrence', 'Metastatic', 'Unknown']


HISTOLOGY_TYPES = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                   'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

def main(options):

    clinical_df = import_metabric_clinical(options.in_clinical)

    clinical_df = filter_clinical_df(clinical_df, options.details_mode)

    # Add Center information
    clinical_df['Center'] = 'METABRIC'

    # Keep only the desired output columns in the desired order
    if options.details_mode:
        clinical_df = clinical_df[OUT_COLS + DETAILS_OUT_COLS]
    else:
        clinical_df = clinical_df[OUT_COLS]

    # Ensure that there is no missing data in any of the columns.
    assert not clinical_df.isnull().values.any()

    # Ensure that there are no duplicated 'Sample_ID' values
    assert not clinical_df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Output the Clinical Data
    clinical_df.to_csv(options.out_clinical, sep="\t", na_rep='', index=False)

    print "##", "-" * 50
    print "## Outfile Summary"
    print "##   Num Samples Clinical Total: %s" % len(clinical_df)
    print "##", "-" * 50

def filter_clinical_df(clinical_df, details_mode):
    """Filter-out clinical samples according to the following logic:
         (1) Exclude Breast_Sarcoma samples (3 in METABRIC dataset)
                - Note: These samples area actually already excluded by condition #1


    Note
    ----
    If 'details_mode:True', then retain these samples but instead just mark the
    'Exclude_Sample' field as 'True' and give reason(s) why the sample should be
    excluded in 'Exclude_Reason' field.
    """

    clinical_df['Exclude_Sample'] = "False"
    clinical_df['Exclude_Reason'] = ""

    # (1) Mark Breast_Sarcoma samples for potential exclusion.
    indices = clinical_df.index[clinical_df.loc[:,'Is_Breast_Sarcoma']]
    mark_row_as_exclude(clinical_df, indices, 'Is_Breast_Sarcoma')

    # (2) If not 'Exclude_Sample', then set 'Exclude_Reason:No_Reason'
    clinical_df.loc[(clinical_df.Exclude_Sample=="False"), 'Exclude_Reason'] = 'No_Reason'

    # Filter-out the samples
    filtered_clinical_df = clinical_df[clinical_df['Exclude_Sample'] == 'False']

    # Consistency checks
    assert filtered_clinical_df['Is_Breast_Carcinoma'].all() == True
    assert filtered_clinical_df['Is_Breast_Sarcoma'].any() == False

    if not details_mode:
        clinical_df = filtered_clinical_df

    return clinical_df

def mark_row_as_exclude(df, indices, new_exclude_reason):
    """Mark the rows in the df that match the specified input indices for
    exclusion and update the 'Exclude_Reason' field."""

    for index in indices:
        so_far_exclude_reason = df.get_value(index, 'Exclude_Reason')

        if so_far_exclude_reason == "":
            so_far_exclude_reason = new_exclude_reason
        else:
            so_far_exclude_reason += "; " + new_exclude_reason

        df.set_value(index, 'Exclude_Sample', 'True')
        df.set_value(index, 'Exclude_Reason', so_far_exclude_reason)

    return df

def import_metabric_clinical(infile):
    """Import the METABRIC clinical file (obtain from cbioportal).

    Notes
    -----
    - Extract the following data:
          (1)  Tumor_Sample_Barcode (e.g. 'MB-7115')
          (2)  ER_Status ('Positive', 'Negative', 'Unknown')
          (3)  PR_Status ('Positive', 'Negative', 'Unknown')
          (4)  HER2_Status ('Positive', 'Negative', 'Unknown')
          (5)  Biopsy_Site_Type ('Primary', 'Local_Recurrence', 'Metastatic', 'Unknown')
          (6)  Histology_Type (see HISTOLOGY_TYPES list)
          (7)  Is_Breast_Sarcoma (True, False)                   [Note: Only for downstream filtering purposes]
          (8)  Is_Breast_Carcinoma (True, False)                 [Note: Only for downstream filtering purposes]

    - Filter-out 'Breast Sarcoma' samples (only 3 counts in METABRIC dataset)
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    assert set(df['Cancer Type'].unique()) == set(['Breast Sarcoma', 'Breast Cancer'])

    df['Is_Breast_Sarcoma'] = (df['Cancer Type'] == 'Breast Sarcoma')  # Not consider a 'Breast Cancer'
    df['Is_Breast_Carcinoma'] = (df['Cancer Type'] == 'Breast Cancer')

    # Check that Patient ID and Sample ID are always the same.
    assert df['Patient ID'].equals(df['Sample ID'])

    # Check that all the require columns exist
    required_columns = ['Sample ID', 'ER Status', 'PR Status', 'HER2 Status', 'Sample Type',
                        'Cancer Type Detailed', 'Is_Breast_Carcinoma', 'Is_Breast_Sarcoma']

    df = df[required_columns]
    df.rename(columns={'Sample ID': 'Tumor_Sample_Barcode'}, inplace=True)
    df.rename(columns={'ER Status': 'ER_Status'}, inplace=True)
    df.rename(columns={'PR Status': 'PR_Status'}, inplace=True)
    df.rename(columns={'HER2 Status': 'HER2_Status'}, inplace=True)
    df.rename(columns={'Sample Type': 'Biopsy_Site_Type'}, inplace=True)
    df.rename(columns={'Cancer Type Detailed': 'Histology_Type'}, inplace=True)

    # Fill missing ER/PR/HER2 status values
    for colname in ['ER_Status', 'PR_Status', 'HER2_Status']:
        df[[colname]] = df[[colname]].fillna(value='Unknown')

    # Map '+' to 'Positive'; map '-' to 'Negative'.
    for colname in ['ER_Status', 'PR_Status', 'HER2_Status']:
        df[colname].replace('+', 'Positive', inplace=True)
        df[colname].replace('-', 'Negative', inplace=True)

        # Ensure that there is no unexpected receptor_status values.
        assert set(df[colname].unique()) == set(['Positive', 'Negative', 'Unknown'])

    # Standardize the values in the 'Histology_Type' column
    df['Histology_Type'].replace('Breast Invasive Ductal Carcinoma',          'Invasive_Ductal_Carcinoma',          inplace=True)
    df['Histology_Type'].replace('Breast Invasive Lobular Carcinoma',         'Invasive_Lobular_Carcinoma',         inplace=True)
    df['Histology_Type'].replace('Breast Mixed Ductal and Lobular Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)
    df['Histology_Type'].replace('Invasive Breast Carcinoma',                 'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Breast Ductal Carcinoma In Situ',           'Other_Breast_Cancer',                inplace=True)
    df['Histology_Type'].replace('Phyllodes Tumor of the Breast',             'Other_Breast_Cancer',                inplace=True)

    df.loc[~df.Is_Breast_Carcinoma, 'Histology_Type'] = 'Not_Breast_Carcinoma'

    # Ensure that there is no unexpected 'Biopsy_Site_Type' values.
    assert set(df['Biopsy_Site_Type'].unique()) <= set(BIOPSY_SITE_TYPES)

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df[df['Is_Breast_Carcinoma']]['Histology_Type'].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    return df

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input METABRIC Clinical Data.")

    parser.add_argument("--details_mode", action="store_true", default=False,
                        help="Output extra rows/columns. For behavior, see header documentation.")

    parser.add_argument("--out_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to output METABRIC Clinical Data.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_clinical: ", repr(options.in_clinical)
    print "##   details_mode:", repr(options.details_mode)
    print "##   out_clinical: ", repr(options.out_clinical)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
