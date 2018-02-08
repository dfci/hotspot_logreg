#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the TCGA BRCA Clinical data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Data source:

        - Supplementary Table 1 of the Cell 2015 Paper:

              http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html

        - Corresponding Cell 2015 data deposited on cBioPortal:

              http://www.cbioportal.org/study?id=brca_tcga_pub2015

    (2) Note on consistency between the 'Exome' column (Tumor_Sample_Barcode) of
        the SI_Table and 'SAMPLE_ID' column of cBioPortal Sample-level Clinical
        data file:

            - Checked that 'Exome' value always start with '{Patient_ID}-01'
              except for TCGA-BH-A1ES --> TCGA-BH-A1ES-06A-12D-A243-09

            - Checked that 'SAMPLE_ID' value always equal'{Patient_ID}-01'.

            - In the discrepancy case ('TCGA-BH-A1ES-06A-12D-A243-09'), manually
              inspected and confirmed that the Histology_Type and Receptor Status
              values of this sample reported in the two sources agrees:
                    Histology_Type: IDC
                    ER_Status:      Positive
                    PR_Status:      Positive
                    HER2_Status:    Negative

EXAMPLES

    ./prep_tcga_clinical.py \
        --in_cbioportal_patient_2015 ../public/tcga_brca/raw/brca_tcga_pub2015/brca_tcga_pub2015/data_bcr_clinical_data_patient.txt \
        --in_cbioportal_sample_2015 ../public/tcga_brca/raw/brca_tcga_pub2015/brca_tcga_pub2015/data_bcr_clinical_data_sample.txt \
        --in_si_table_cell_2015 ../public/clinical/raw/tcga_brca/Ciriello_Cell_2015.Table_S1.xlsx \
        --out_clinical ../public/clinical/tcga.clinical_data.tsv

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (April-2017)
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

HISTOLOGY_TYPES = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                   'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

GENDER_TYPES = ['Female', 'Male', 'Unknown']

def main(options):

    cli_patient_df = import_cbioportal_patient_data(options.in_cbioportal_patient_2015)

    cli_sample_df = import_cbioportal_sample_data(options.in_cbioportal_sample_2015)

    cli_si_table_df = import_si_table_cell_2015_cli_data(options.in_si_table_cell_2015)

    clinical_df = merge_clinical_df(cli_patient_df, cli_sample_df, cli_si_table_df)

    # Assume that all TCGA-BRCA samples have 'Biopsy_Site_Type:Primary'
    clinical_df['Biopsy_Site_Type'] = 'Primary'

    # Add Center information
    clinical_df['Center'] = 'TCGA'

    # Keep only the desired output columns in the desired order.
    clinical_df = clinical_df[OUT_COLS]

    # Ensure that there is no missing data in any of the columns.
    assert not clinical_df.isnull().values.any()

    # Output the Clinical Data
    clinical_df.to_csv(options.out_clinical, sep="\t", na_rep='', index=False)

    print "##", "-" * 50
    print "## Outfile Summary"
    print "##   Total # Samples (Clinical): %s" % len(clinical_df)
    print "##", "-" * 50

def merge_clinical_df(cli_patient_df, cli_sample_df, cli_si_table_df):
    """Perform outer join on the three input dfs using 'Patient_ID' as the join
    key.


    Note
    ----
    (1) Will ensure that the set of 'Patient_ID' found in the three input DFs
        are the same. So will obtain same result if switch to using INNER JOIN.

    """

    # Check that there is no overlap columns between the three input dfs (except
    # for the 'Patient_ID' join key)
    colnames_so_far = set()
    for colname in list(cli_patient_df) + list(cli_sample_df) + list(cli_si_table_df):
        if colname != 'Patient_ID':
            if colname in colnames_so_far:
                raise Exception("Found duplicated colname (%s)." % colname)
            colnames_so_far.add(colname)

    assert set(cli_patient_df['Patient_ID'].unique()) == set(cli_sample_df['Patient_ID'].unique())
    assert set(cli_patient_df['Patient_ID'].unique()) == set(cli_si_table_df['Patient_ID'].unique())

    # Perform outer joins
    tmp_df = pd.merge(cli_patient_df, cli_sample_df,
                      how="outer",  # OUTER JOIN
                      on="Patient_ID",
                      sort=False,
                      indicator='indicator_column1')

    df = pd.merge(tmp_df, cli_si_table_df,
                  how="outer",  # OUTER JOIN
                  on="Patient_ID",
                  sort=False,
                  indicator='indicator_column2')

    # Ensure that is is no 'orphan' Patient_ID that is found in only some of the
    # input DFs.
    assert set(df['indicator_column1'].unique()) == set(['both'])
    assert set(df['indicator_column2'].unique()) == set(['both'])

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Consistency checks
    # Ensure that there is no unexpected receptor_status values.
    for colname in RECEPTOR_STATUS_COLS:
        assert set(df[colname].unique()) == set(RECEPTOR_STATUSES)

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df['Histology_Type'].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there is no unexpected 'Gender' values.
    assert set(df['Gender'].unique()) <= set(GENDER_TYPES)

    # Ensure that there are no duplicated 'Patient_ID' values
    assert not df.duplicated(subset=['Patient_ID']).any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    return df

def import_si_table_cell_2015_cli_data(infile):
    """Import the TCGA-BRCA clinical data from Supplemental Table 1 of the 2015
    Cell paper

    Notes
    -----
    - Extract the following data in standardized format:
        (1)  Patient_ID (e.g. 'TCGA-A2-A0T2')
        (2)  Tumor_Sample_Barcode (e.g. 'TCGA-A2-A0T2-01A-11W-A097-09')
        (3)  Histology_Type (see HISTOLOGY_TYPES list)
        (4)  ER_Status_SI_Table_NOT_USED (see RECEPTOR_STATUSES list)
        (5)  PR_Status_SI_Table_NOT_USED (see RECEPTOR_STATUSES list)
        (6)  HER2_Status_SI_Table_NOT_USED (see RECEPTOR_STATUSES list)
    """

    df = pd.read_excel(infile, sheetname="Suppl. Table 1", skiprows=2)

    # Select all the require columns + implicitly check that all required columns exist
    required_columns = ['Case.ID', 'Exome', 'Final Pathology', 'ER IHC', 'PR IHC', 'HER2 IHC']
    df = df[required_columns]
    df.rename(columns={'Case.ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'Exome': 'Tumor_Sample_Barcode'}, inplace=True)
    df.rename(columns={'Final Pathology': 'Histology_Type'}, inplace=True)
    df.rename(columns={'ER IHC': 'ER_Status_SI_Table_NOT_USED'}, inplace=True)
    df.rename(columns={'PR IHC': 'PR_Status_SI_Table_NOT_USED'}, inplace=True)
    df.rename(columns={'HER2 IHC': 'HER2_Status_SI_Table_NOT_USED'}, inplace=True)

    # Standardize the values in the 'Histology_Type' column
    df['Histology_Type'].replace('IDC',           'Invasive_Ductal_Carcinoma',          inplace=True)
    df['Histology_Type'].replace('ILC',           'Invasive_Lobular_Carcinoma',         inplace=True)
    df['Histology_Type'].replace('Mixed.IDC.ILC', 'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)
    df['Histology_Type'].replace('Other',         'Other_Invasive_Breast_Carcinoma',    inplace=True)

    # Fill missing ER/PR/HER2 status values
    for colname in RECEPTOR_STATUS_COLS:
        colname = colname + '_SI_Table_NOT_USED'
        df[[colname]] = df[[colname]].fillna(value='Unknown')

    # Map of possible missing/equivocal results to 'Unknown'
    for colname in RECEPTOR_STATUS_COLS:
        colname = colname + '_SI_Table_NOT_USED'
        df[colname].replace('[Not Evaluated]', 'Unknown', inplace=True)
        df[colname].replace('[Not Available]', 'Unknown', inplace=True)
        df[colname].replace('Equivocal', 'Unknown', inplace=True)
        df[colname].replace('Indeterminate', 'Unknown', inplace=True)
        df[colname].replace('#N/A', 'Unknown', inplace=True)

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Patient_ID' values
    assert not df.duplicated(subset=['Patient_ID']).any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Check that for all rows, Tumor_Sample_Barcode has Patient_ID + '-01' as its prefix.
    for index, row in df.iterrows():
        if row['Patient_ID'] != 'TCGA-BH-A1ES':
            # See header documentation for rationale of this check.
            assert row['Tumor_Sample_Barcode'].startswith(row['Patient_ID'] + '-01')
        else:
            assert row['Tumor_Sample_Barcode'].startswith(row['Patient_ID'] + '-06')

    return df

def import_cbioportal_sample_data(infile):
    """Import the Sample-level Clinical Data for the 2015 Cell TCGA dataset
    (downloaded from cBioPortal).

    Notes
    -----
    - Extract the following data in standardized format:
        (1)  Patient_ID (e.g. 'TCGA-LQ-A4E4')
        (2)  ER_Status (see RECEPTOR_STATUSES list)
        (3)  PR_Status (see RECEPTOR_STATUSES list)
        (4)  HER2_Status (see RECEPTOR_STATUSES list)
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    remove_duplicate_rows(df)

    # Check that for all rows, 'SAMPLE_ID is always equal to PATIENT_ID + '-01'
    for index, row in df.iterrows():
        assert row['SAMPLE_ID'].startswith(row['PATIENT_ID'] + '-01')

    # Select all the require columns + implicitly check that all required columns exist
    required_columns = ['PATIENT_ID', 'SAMPLE_ID', 'ER_STATUS_BY_IHC', 'PR_STATUS_BY_IHC', 'IHC_HER2', 'HER2_FISH_STATUS']
    df = df[required_columns]
    df.rename(columns={'PATIENT_ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'SAMPLE_ID': 'Sample_ID'}, inplace=True)
    df.rename(columns={'ER_STATUS_BY_IHC': 'ER_Status'}, inplace=True)
    df.rename(columns={'PR_STATUS_BY_IHC': 'PR_Status'}, inplace=True)
    df.rename(columns={'IHC_HER2': 'HER2_IHC'}, inplace=True)
    df.rename(columns={'HER2_FISH_STATUS': 'HER2_FISH'}, inplace=True)

    df = infer_HER2_status(df)
    df = df.drop(['HER2_IHC', 'HER2_FISH'], 1)

    # Map of possible missing/equivocal results to 'Unknown'
    for colname in RECEPTOR_STATUS_COLS:
        df[colname].replace('[Not Available]', 'Unknown', inplace=True)
        df[colname].replace('Indeterminate', 'Unknown', inplace=True)

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Patient_ID' values
    assert not df.duplicated(subset=['Patient_ID']).any()

    return df

def remove_duplicate_rows(df):
    """There appear to be 13 'completely' duplicate rows in the cBioPortal
    sample-level clinical data file. Remove these duplicate rows from
    the dataframe."""

    exclude_indices_keep_false = df.index[df.duplicated(keep=False)]

    exclude_indices_keep_first = df.index[df.duplicated(keep='first')]

    assert len(exclude_indices_keep_false) == 26
    assert len(exclude_indices_keep_first) == 13

    print "##"
    print "## WARNING: Keep only the first instance of the following duplicated",
    print "rows from Sample-level Clinical Data File (cBioPortal; Cell 2015 TCGA-BRCA):"""
    print "##"
    print df[df.index.isin(exclude_indices_keep_false)][['PATIENT_ID', 'SAMPLE_ID', 'OTHER_SAMPLE_ID']]

    df.drop(exclude_indices_keep_first, inplace=True)

    return df

def infer_HER2_status(df):
    """Infer the sample's HER2_Status from 'HER2_IHC' and 'HER2_FISH' data

    Notes
    -----

    (1) Here is counts of 'HER2_IHC' and 'HER2_FISH' for the Cell 2015 dataset
        found in the Sample-level Clinical data downloaded from cBioPortal
        (April 05th, 2017).

        Command: df.groupby(['HER2_IHC','HER2_FISH']).size()

        (A) Post-Standardize Value:

            HER2_IHC   HER2_FISH    Counts         Inferred_HER2_Status
            -----------------------------------------------------------
            Equivocal  Equivocal         3   -->   Unknown
                       Negative        109   -->   Negative
                       Positive         16   -->   Positive
                       Unknown          12   -->   Unknown
            Negative   Negative         93   -->   Negative
                       Positive          2   -->   Positive     [Note: Inconsistent!]
                       Unknown         330   -->   Negative
            Positive   Equivocal         2   -->   Positive
                       Negative          7   -->   Positive     [Note: Inconsistent!]
                       Positive         34   -->   Positive
                       Unknown          78   -->   Positive
            Unknown    Negative         42   -->   Negative
                       Positive          9   -->   Positive
                       Unknown          93   -->   Unknown
    """

    for colname in ['HER2_IHC', 'HER2_FISH']:
        df[colname].replace('[Not Available]', 'Unknown',  inplace=True)
        df[colname].replace('Indeterminate', 'Unknown',  inplace=True)
        assert set(df[colname].unique()) == set(['Positive', 'Negative', 'Equivocal', 'Unknown'])

    df['HER2_Status'] = "Unknown"
    df.loc[(df.HER2_IHC=='Negative'), 'HER2_Status'] = 'Negative'
    df.loc[(df.HER2_FISH=='Negative'), 'HER2_Status'] = 'Negative'
    df.loc[(df.HER2_IHC=='Positive'), 'HER2_Status'] = 'Positive'  # Overide the Negative
    df.loc[(df.HER2_FISH=='Positive'), 'HER2_Status'] = 'Positive'  # Overide the Negative

    return df

def import_cbioportal_patient_data(infile):
    """Import the Patient-level Clinical Data for the 2015 Cell TCGA dataset
    (downloaded from cBioPortal).

    Notes
    -----
    - Extract the following data in standardized format:
        (1)  Patient_ID (e.g. 'TCGA-A2-A0T2')
        (2)  Histology_Type_cBioPortal_NOT_USED (see HISTOLOGY_TYPES list)
        (3)  Gender (see GENDER_TYPES list)
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Select all the require columns + implicitly check that all required columns exist
    required_columns = ['PATIENT_ID', 'HISTOLOGICAL_DIAGNOSIS', 'GENDER']
    df = df[required_columns]

    # Standardize the column names
    df.rename(columns={'PATIENT_ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'HISTOLOGICAL_DIAGNOSIS': 'Histology_Type_cBioPortal_NOT_USED'}, inplace=True)
    df.rename(columns={'GENDER': 'Gender'}, inplace=True)

    colname = 'Histology_Type_cBioPortal_NOT_USED'

    # Standardize the values in the 'Histology_Type' column
    df[colname].replace('Infiltrating Ductal Carcinoma',    'Invasive_Ductal_Carcinoma',          inplace=True)
    df[colname].replace('Infiltrating Lobular Carcinoma',   'Invasive_Lobular_Carcinoma',         inplace=True)
    df[colname].replace('Mixed Histology (please specify)', 'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)
    df[colname].replace('Infiltrating Carcinoma NOS',       'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[colname].replace('Medullary Carcinoma',              'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[colname].replace('Mucinous Carcinoma',               'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[colname].replace('Metaplastic Carcinoma',            'Other_Breast_Cancer',                inplace=True)
    df[colname].replace('Other, specify',                   'Unknown',                            inplace=True)
    df[colname].replace('[Not Available]',                  'Unknown',                            inplace=True)

    # Standardize the values in the 'Gender' column
    df['Gender'] = df['Gender'].fillna(value='Unknown')
    df['Gender'].replace('FEMALE', 'Female', inplace=True)
    df['Gender'].replace('MALE', 'Male', inplace=True)

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Patient_ID' values
    assert not df.duplicated(subset=['Patient_ID']).any()

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

    parser.add_argument("--in_cbioportal_patient_2015", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Patient-Level Clinical Data of the Cell 2015 paper (from cBioportal).")

    parser.add_argument("--in_cbioportal_sample_2015", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Sample-Level Clinical Data of the Cell 2015 paper (from cBioportal).")

    parser.add_argument("--in_si_table_cell_2015", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Clinical Data of the Cell 2015 paper (from the SI Table).")

    parser.add_argument("--out_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to output Clinical Data.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_cbioportal_patient_2015: ", repr(options.in_cbioportal_patient_2015)
    print "##   in_cbioportal_sample_2015: ", repr(options.in_cbioportal_sample_2015)
    print "##   in_si_table_cell_2015: ", repr(options.in_si_table_cell_2015)
    print "##   out_clinical: ", repr(options.out_clinical)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
