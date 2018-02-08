#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the DFCI Breast Cancer Clinical data into a format that is
    suitable for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Support 'details' mode, which will:

          (A) Include 'Panel_Version' column

          (B) Retain samples which will be normally filtered out (e.g. duplicated
              samples from a patient, Breast Sarcoma samples) and instead
              add 'Exclude_Sample' column  ('True', 'False') and 'Exclude_Reason'
              column (text describing the reason to exclude sample).

              For sample exclusion logic, see filter_clinical_df() function.

        Note
        ----
        Purpose of the 'details' mode is to include extra rows/columns in the
        preprocessed clinical outfile that can be consumed by the 'prep_dfci_maf.py'
        script BUT is not used/supported by the actual mutation counts pipeline.

EXAMPLES

    # (1) Generate the 'regular' preprocessed DFCI clinical outfile.
    ./prep_dfci_clinical.py \
        --in_camd_clinical ../private/raw/06Feb2017/data_samples.txt \
        --in_camd_gender ../private/raw/06Feb2017/data_patients.txt \
        --in_haoguo_clinical ../private/raw/11Apr2017/clinical_de_identified.sample.txt \
        --out_clinical ../private/clinical/dfci.clinical_data.r5.check.tsv

    # (2) Generate the 'details' preprocessed DFCI clinical outfile.
    ./prep_dfci_clinical.py \
        --in_camd_clinical ../private/raw/06Feb2017/data_samples.txt \
        --in_camd_gender ../private/raw/06Feb2017/data_patients.txt \
        --in_haoguo_clinical ../private/raw/11Apr2017/clinical_de_identified.sample.txt \
        --details_mode \
        --out_clinical ../private/clinical/dfci.clinical_data.r5.check.details.tsv


AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (April-2017)
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

OUT_COLS = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status', 'HER2_Status', 'Biopsy_Site_Type', 'Histology_Type', 'Gender']
DETAILS_OUT_COLS = ['Panel_Version', 'Exclude_Sample', 'Exclude_Reason']

RECEPTOR_STATUS_COLS = ['ER_Status', 'PR_Status', 'HER2_Status']
RECEPTOR_STATUSES = ['Positive', 'Negative', 'Unknown']

BIOPSY_SITE_TYPES = ['Primary', 'Local_Recurrence', 'Metastatic', 'Unknown']

HISTOLOGY_TYPES = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                   'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

GENDER_TYPES = ['Female', 'Male', 'Unknown']

PANEL_VERSIONS = ['1', '2', '3']

def main(options):

    camd_sample_df = import_camd_clinical(options.in_camd_clinical)

    camd_gender_df = import_camd_gender(options.in_camd_gender)

    camd_clinical_df = merge_camd_df(camd_sample_df, camd_gender_df)

    haoguo_clinical_df = import_haoguo_clinical(options.in_haoguo_clinical)

    clinical_df = merge_clinical_df(camd_clinical_df, haoguo_clinical_df)

    check_merged_clinical_df(clinical_df)

    # For debugging/logging purposes
    print_patients_with_inconsistent_cancer_types(clinical_df)

    clinical_df = filter_clinical_df(clinical_df, options.details_mode)

    # Set Center to 'DFCI'
    clinical_df['Center'] = 'DFCI'

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
         (1) Include only Is_Breast_Carcinoma (Exclude all other Cancer Types)

         (2) Exclude Breast_Sarcoma samples (3 in DFCI dataset)
                - Note: These samples area actually already excluded by condition #1

         (3) Select and retain only 1 sample per patient

    Note
    ----
    If 'details_mode:True', then retain these samples but instead just mark the
    'Exclude_Sample' field as 'True' and give reason(s) why the sample should be
    excluded in 'Exclude_Reason' field.
    """

    clinical_df['Exclude_Sample'] = "False"
    clinical_df['Exclude_Reason'] = ""

    # (1) Mark non-'Breast_Carcinoma' samples for potential exclusion.
    indices = clinical_df.index[clinical_df.loc[:,'Is_Breast_Carcinoma'] == False]
    mark_row_as_exclude(clinical_df, indices, 'Not_Breast_Carcinoma')

    # (2) Mark Breast_Sarcoma samples for potential exclusion.
    indices = clinical_df.index[clinical_df.loc[:,'Is_Breast_Sarcoma']]
    mark_row_as_exclude(clinical_df, indices, 'Is_Breast_Sarcoma')

    # (3) Mark extra samples of a patient for potential exclusion.
    clinical_df = mark_extra_samples_of_patient(clinical_df)

    # (4) If not 'Exclude_Sample', then set 'Exclude_Reason:No_Reason'
    clinical_df.loc[(clinical_df.Exclude_Sample=="False"), 'Exclude_Reason'] = 'No_Reason'

    # Filter-out the samples
    filtered_clinical_df = clinical_df[clinical_df['Exclude_Sample'] == 'False']

    # Consistency checks
    assert filtered_clinical_df['Is_Breast_Carcinoma'].all() == True
    assert filtered_clinical_df['Is_Breast_Sarcoma'].any() == False
    assert (filtered_clinical_df["Patient_ID"].value_counts() > 1).any() == False

    if not details_mode:
        clinical_df = filtered_clinical_df

    return clinical_df

def mark_extra_samples_of_patient(clinical_df):
    """If a patient has multiple breast cancer samples sequenced, then mark
    the 'extra' samples for potential downstream exclusion.

    Choose the 'retained' samples that is 'unmarked' using the following citeria:

        - Fistly, choose sample that had not already been marked for exclusion
          for other reasons (e.g. not a 'Breast Cancer sample')

        - Secondly, choose 'Biopsy_Site_Type' in the following order:
                (1) 'Primary'              [First]
                (2) 'Local_Recurrence'
                (3) 'Metastatic'
                (4) 'Unknown'              [Last]

        - Thirdly, choose the sample with the earlier 'Age_At_Seq_Report.

        - Fourthly, if there is still a tie, then make the choice deterministic by
          choosing by alphanumeric sorting order (e.g. choose
          'CBIO_P1008_S1' over 'CBIO_P1008_S2').
    """

    # Create dictionary mapping from patient to sample + associated information.
    patient2samples = dict()
    for index, row in clinical_df.iterrows():
        patient_id = row["Patient_ID"]
        sample_id = row["Tumor_Sample_Barcode"]
        biopsy_site_type = row['Biopsy_Site_Type']

        age = float(row['Age_At_Seq_Report'])
        age_times_100 = age * 100.0
        assert (age_times_100).is_integer
        age_times_100 = int(age_times_100)

        already_marked = row["Exclude_Sample"]

        base_multipler = 10000

        # Consistency checks
        assert biopsy_site_type in BIOPSY_SITE_TYPES
        assert already_marked in ['True', 'False']
        assert age > 0.0
        assert age < 100.0
        assert age_times_100 > 1
        assert age_times_100 < base_multipler

        # Compute sample_priority (where lower value equals higher priority)
        sample_priority = 0
        if already_marked == 'False':
            sample_priority += 0 * 10**3 * base_multipler
        elif already_marked == 'True':
            sample_priority += 1 * 10**3 * base_multipler
        else:
            raise Exception("Invalid 'Already_Marked' value.")

        if biopsy_site_type == 'Primary':
            sample_priority += 0 * 10**1 * base_multipler
        elif biopsy_site_type == 'Local_Recurrence':
            sample_priority += 1 * 10**1 * base_multipler
        elif biopsy_site_type == 'Metastatic':
            sample_priority += 2 * 10**1 * base_multipler
        elif biopsy_site_type == 'Unknown':
            sample_priority += 3 * 10**1 * base_multipler
        else:
            raise Exception("Invalid 'Biopsy_Site_Type' type.")

        sample_priority += 1 * 10**0 * age_times_100

        assert isinstance(sample_priority, int)

        if patient_id not in patient2samples:
            patient2samples[patient_id] = dict()
        assert sample_id not in patient2samples[patient_id]
        patient2samples[patient_id][sample_id] = sample_priority

    # Create a list of samples to keep.
    keep_samples = set()
    for patient in patient2samples:
        samples_dict = patient2samples[patient]

        # Get Sample(s) with the minimum priority
        min_priority = None
        for sample in samples_dict:
            sample_priority = samples_dict[sample]
            if min_priority == None:
                min_priority = sample_priority

            min_priority = min(sample_priority, min_priority)

        samples = [sample for sample in samples_dict.keys() if samples_dict[sample] == min_priority]
        keep_sample = sorted(samples)[0]  # Alphanumeric Sorting Order for tie breaking
        assert sample not in keep_samples
        keep_samples.add(keep_sample)

    # Mark the extra samples of a patient for potential exclusion.
    indices = clinical_df.index[clinical_df['Tumor_Sample_Barcode'].isin(keep_samples) == False]
    mark_row_as_exclude(clinical_df, indices, 'Keep_Only_1_Sample_Per_Patient')

    # Print the Marked Breast_Carcinoma samples to stdout for debugging/logging purposes
    marked_brca_df = clinical_df[(~clinical_df['Tumor_Sample_Barcode'].isin(keep_samples)) & (clinical_df['Is_Breast_Carcinoma'])]
    print "##"
    print "## DEBUG: Marked %s Breast_Carcinoma Samples for exclusion" % len(marked_brca_df),
    print "from %s multisamples Patients:" % marked_brca_df['Patient_ID'].nunique()
    print
    print marked_brca_df[["Patient_ID", "Tumor_Sample_Barcode", "Biopsy_Site_Type", "Is_Breast_Carcinoma", "Histology_Type"]]
    print

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

def print_patients_with_inconsistent_cancer_types(clinical_df):
    """For debugging purposes, print a list of patients with:
        (1) Two or more sample sequenced
        (2) Inconsistent Cancer_Type between the patient's samples.
        (3) At least 1 of the sample's Cancer_Type is 'Breast_Carcinoma'
    """

    patients_to_print = set()

    # Create dictionary mapping from patient to sample + associated information.
    patient2samples = dict()
    for index, row in clinical_df.iterrows():
        patient_id = row['Patient_ID']
        sample_id = row['Tumor_Sample_Barcode']
        is_breast_carcinoma = row['Is_Breast_Carcinoma']

        if patient_id not in patient2samples:
            patient2samples[patient_id] = dict()

        assert sample_id not in patient2samples[patient_id]

        patient2samples[patient_id][sample_id] = dict()
        patient2samples[patient_id][sample_id]['Is_Breast_Carcinoma'] = is_breast_carcinoma

    for patient in patient2samples:
        samples_dict = patient2samples[patient]

        has_breast_cancer_sample = False
        is_breast_carcinoma_vals = set()
        for sample_id in samples_dict:
            is_breast_carcinoma = samples_dict[sample_id]['Is_Breast_Carcinoma']
            is_breast_carcinoma_vals.add(is_breast_carcinoma)

        if True in is_breast_carcinoma_vals:
            if len(is_breast_carcinoma_vals) > 1:
                patients_to_print.add(patient)

    tmp_df = clinical_df[clinical_df["Patient_ID"].isin(patients_to_print)]
    print "## DEBUG: Found %s patients with inconsistent Cancer Types:" % len(patients_to_print)
    print
    print tmp_df[["Patient_ID", "Tumor_Sample_Barcode", "Biopsy_Site_Type", "CANCER_TYPE", "Histology_Type"]]
    print

def merge_camd_df(camd_sample_df, camd_gender_df):
    """Merge the sample-level CAMD Clinical DF and the CAMD Gender DF."""

    # Check that there is no overlap columns between the two input CAMD DF
    # (except for the 'Patient_ID' join key)
    non_key_camd_sample = set(camd_sample_df) - set(['Patient_ID'])
    non_key_camd_gender = set(camd_gender_df) - set(['Patient_ID'])
    assert len(non_key_camd_sample & non_key_camd_gender) == 0

    # Perform outer join
    df = pd.merge(camd_sample_df, camd_gender_df,
                  how="outer",  # OUTER JOIN
                  on="Patient_ID",
                  sort=False,
                  indicator='indicator_column1')

    # Ensure that is is no 'orphan' Patient_ID that is found in only 1 of the 2
    # input CAMD DFs.
    assert set(df['indicator_column1'].unique()) <= set(['both'])

    return df

def merge_clinical_df(camd_cli_df, haoguo_cli_df):
    """Merge CAMD Clinical DF and the Hao-Guo Clinical DF."""

    # Ensure consistency between the 'Patient_ID' and 'Tumor_Sample_Barcode'
    # matches.
    check_patient_and_sample_matches(camd_cli_df, haoguo_cli_df)
    haoguo_cli_df = haoguo_cli_df.drop(['Patient_ID'], 1)

    # Check that there is no overlap columns between the two input CAMD DF
    # (except for the 'Tumor_Sample_Barcode' join key)
    non_key_camd_cli = set(camd_cli_df) - set(['Tumor_Sample_Barcode'])
    non_key_haoguo_cli = set(haoguo_cli_df) - set(['Tumor_Sample_Barcode'])
    assert len(non_key_camd_cli & non_key_haoguo_cli) == 0

    # Perform left outer join.
    df = pd.merge(camd_cli_df, haoguo_cli_df,
                  how="left",  # LEFT OUTER JOIN
                  on="Tumor_Sample_Barcode",
                  sort=False,
                  indicator='indicator_column2')

    # This should be True since performing left_outer join
    assert set(df['indicator_column2'].unique()) <= set(['both', 'left_only'])

    # Fill-in missing values in Column derived from the Hao-Guo Clinical dataset
    assert non_key_haoguo_cli == set(RECEPTOR_STATUS_COLS + ['Histology_Type_HAOGUO'])
    for colname in non_key_haoguo_cli:
        df[colname] = df[colname].fillna(value='Unknown')

    df = combine_histology_type_data(df)
    df = df.drop(['Histology_Type_CAMD', 'Histology_Type_HAOGUO'], 1)

    return df

def check_merged_clinical_df(df):
    """Check the contents of the Merged Clinical DF."""

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there is no unexpected receptor_status values.
    for colname in RECEPTOR_STATUS_COLS:
        assert set(df[colname].unique()) == set(RECEPTOR_STATUSES)

    # Ensure that there is no unexpected 'Biopsy_Site_Type' values.
    assert set(df['Biopsy_Site_Type'].unique()) <= set(BIOPSY_SITE_TYPES)

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df[df['Is_Breast_Carcinoma']]['Histology_Type'].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there is no unexpected 'Gender' values.
    assert set(df['Gender'].unique()) <= set(GENDER_TYPES)

    # Ensure that there is no unexpected 'Panel_Version' values.
    assert set(df['Panel_Version'].unique()) <= set(PANEL_VERSIONS)

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    return df

def combine_histology_type_data(cli_df):
    """Create new 'Histology_Type' column which merges the histology_type data
    from the CAMD's and HAO-GUO's data sources."""

    # Ensure that there is no missing data in any of the columns.
    assert not cli_df.isnull().values.any()

    assert set(cli_df[cli_df['Is_Breast_Carcinoma']]['Histology_Type_CAMD'].unique()) <= set(HISTOLOGY_TYPES)
    assert set(cli_df['Histology_Type_HAOGUO'].unique()) <= set(HISTOLOGY_TYPES)

    cli_df['Histology_Type'] = cli_df['Histology_Type_CAMD']

    # If the Hao-Gao's dataset assigns a known 'histology_type' value to the
    # sample, then this assignment take precedence over the CAMD's assignment.
    indices = cli_df.index[cli_df.loc[:,'Histology_Type_HAOGUO'] != 'Unknown']

    for index in indices:
        histology_type = cli_df.get_value(index, 'Histology_Type_HAOGUO')
        cli_df.set_value(index, 'Histology_Type', histology_type)

    return cli_df

def check_patient_and_sample_matches(camd_cli_df, haoguo_cli_df):
    """Ensure consistency between the 'Patient_ID' and 'Tumor_Sample_Barcode'
    matches.
    """

    match_patients = set(camd_cli_df['Patient_ID'].unique()) & set(haoguo_cli_df['Patient_ID'].unique())
    match_samples = set(camd_cli_df['Tumor_Sample_Barcode'].unique()) & set(haoguo_cli_df['Tumor_Sample_Barcode'].unique())

    sample2patient_camd, patient2samples_camd = map_sample_and_patient(camd_cli_df)
    sample2patient_haoguo, patient2samples_haoguo = map_sample_and_patient(haoguo_cli_df)

    # For each 'Tumor_Sample_Barcode' match, there is a corresponding 'Patient_ID'
    # match.
    for sample_id in match_samples:
        patient_camd = sample2patient_camd[sample_id]
        patient_haoguo = sample2patient_haoguo[sample_id]
        assert patient_camd == patient_haoguo

    # For each 'Patient_ID' match, there is a corresponding 'Tumor_Sample_Barcode'
    # match.
    for patient_id in match_patients:
        samples_camd = patient2samples_camd[patient_id]
        sample_haoguo = patient2samples_haoguo[patient_id]
        assert len(samples_camd & sample_haoguo) > 0

    # Same number of matches.
    assert len(match_patients) == len(match_samples)

def map_sample_and_patient(cli_df):
    """Create a dictionary mapping 'Tumor_Sample_Barcode' and 'Patient_ID' """

    sample2patient = dict()
    patient2samples = dict()

    for index, row in cli_df.iterrows():
        sample = row['Tumor_Sample_Barcode']
        patient = row['Patient_ID']

        if sample not in sample2patient:
            sample2patient[sample] = set()
        sample2patient[sample].add(patient)

        if patient not in patient2samples:
            patient2samples[patient] = set()
        patient2samples[patient].add(sample)

    for sample in sample2patient:
        patients = sample2patient[sample]
        assert len(patients) == 1
        sample2patient[sample] = list(patients)[0]

    return sample2patient, patient2samples

def import_camd_clinical(infile):
    """Import the DFCI Sample-level Clinical File [obtained from the Center for
    Advanced Molecular Diagnostics (CAMD) lab].

    Notes
    -----
    - Extract the following data:
          (1)  Patient_ID (e.g. 'CBIO_P10001')                   [Note: Used only for joining clinical tables + logging]
          (2)  Tumor_Sample_Barcode (e.g. 'CBIO_P10001_S1')
          (3)  Biopsy_Site_Type (see BIOPSY_SITE_TYPES list)
          (4)  Histology_Type_CAMD (see HISTOLOGY_TYPES list)
          (5)  Panel_Version (see PANEL_VERSIONS list)           [Note: Only output in details_mode]
          (6)  CANCER_TYPE                                       [Note: Only for debugging/logging purposes]
          (6)  Is_Breast_Sarcoma (True, False)                   [Note: Only for downstream filtering purposes]
          (7)  Is_Breast_Carcinoma (True, False)                 [Note: Only for downstream filtering purposes]
          (8)  Age_At_Seq_Report                                 [Note: Only for downstream filtering purposes]
    """
    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    df['Is_Breast_Sarcoma'] = (df['CANCER_TYPE'] == 'Breast Sarcoma')  # Not consider a 'Breast Cancer'
    df['Is_Breast_Carcinoma'] = (df['CANCER_TYPE'] == 'Breast Carcinoma')

    # Retain only the required columns
    required_columns = ['PATIENT_ID', 'SAMPLE_ID', 'ONCOTREE_BIOPSY_SITE_TYPE', 'CANCER_TYPE_DETAILED', 'PANEL_VERSION', 'CANCER_TYPE',
                        'Is_Breast_Carcinoma', 'Is_Breast_Sarcoma', 'AGE']

    assert set(required_columns) <= set(df)
    df = df[required_columns]

    histology_colname = 'Histology_Type_CAMD'

    df.rename(columns={'PATIENT_ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'SAMPLE_ID': 'Tumor_Sample_Barcode'}, inplace=True)
    df.rename(columns={'ONCOTREE_BIOPSY_SITE_TYPE': 'Biopsy_Site_Type'}, inplace=True)
    df.rename(columns={'CANCER_TYPE_DETAILED': histology_colname}, inplace=True)
    df.rename(columns={'PANEL_VERSION': 'Panel_Version'}, inplace=True)
    df.rename(columns={'AGE': 'Age_At_Seq_Report'}, inplace=True)

    # Standardize the values in the 'Biopsy_Site_type' column
    df['Biopsy_Site_Type'].replace('Local Recurrence', 'Local_Recurrence', inplace=True)
    df['Biopsy_Site_Type'].replace('Metastatic Recurrence', 'Metastatic', inplace=True)
    df['Biopsy_Site_Type'].replace('Any/Other', 'Unknown', inplace=True)
    df['Biopsy_Site_Type'].replace('Unspecified', 'Unknown', inplace=True)
    df['Biopsy_Site_Type'].replace('Not Applicable', 'Unknown', inplace=True)

    # Standardize the values in the histology_type column
    df[histology_colname].replace('Breast Invasive Ductal Carcinoma',          'Invasive_Ductal_Carcinoma',          inplace=True)
    df[histology_colname].replace('Breast Invasive Lobular Carcinoma',         'Invasive_Lobular_Carcinoma',         inplace=True)
    df[histology_colname].replace('Breast Mixed Ductal and Lobular Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)
    df[histology_colname].replace('Invasive Breast Carcinoma',                 'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Adenoid Cystic Breast Cancer',              'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Breast Invasive Carcinosarcoma, NOS',       'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Breast Invasive Cancer, NOS',               'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Breast Carcinoma with Signet Ring',         'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Breast Invasive Mixed Mucinous Carcinoma',  'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Solid Papillary Carcinoma of the Breast',   'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df[histology_colname].replace('Adenomyoepithelioma of the Breast',         'Other_Breast_Cancer',                inplace=True)
    df[histology_colname].replace('Breast Ductal Carcinoma In Situ',           'Other_Breast_Cancer',                inplace=True)
    df[histology_colname].replace('Carcinoma with Chondroid Metaplasia',       'Other_Breast_Cancer',                inplace=True)
    df[histology_colname].replace('Metaplastic Breast Cancer',                 'Other_Breast_Cancer',                inplace=True)
    df[histology_colname].replace('Paget Disease of the Nipple',               'Other_Breast_Cancer',                inplace=True)

    df.loc[~df.Is_Breast_Carcinoma, histology_colname] = 'Not_Breast_Carcinoma'

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Check that for all rows, Exome_ID has Sample_ID as its prefix.
    for index, row in df.iterrows():
        assert row['Tumor_Sample_Barcode'].startswith(row['Patient_ID'] + '_S')

    return df

def import_camd_gender(infile):
    """Import the DFCI Patients' gender data file [obtained from the Center for
    Advanced Molecular Diagnostics (CAMD) lab].

    Notes
    -----
    - Extract the following data:
          (1)  PATIENT_ID (e.g. 'CBIO_P10001')  [Note: Used only for joining clinical tables]
          (2)  Gender (see GENDER_TYPES list)
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Retain only the required columns
    required_columns = ['PATIENT_ID', 'GENDER']

    df.rename(columns={'PATIENT_ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'GENDER': 'Gender'}, inplace=True)

    # Ensure that there are no duplicated 'Patient_ID' values
    assert not df.duplicated(subset=['Patient_ID']).any()

    return df

def import_haoguo_clinical(infile):
    """Import the DFCI clinical data file obtained from Hao Guo (BCB).

    For meaning of values in each column, also see the accompanying dictionary
    document (located in the same directory as the Clinical Data file):

        Clinicaldata_forIntelCCCproject_4-10-2017_Dictionary.xls

    Notes
    -----
    - Extract the following data:
          (1)  Patient_ID (e.g. 'CBIO_P10001')                           [Note: Used only for consistency_check]
          (2)  Tumor_Sample_Barcode (e.g. 'CBIO_P10001_S1')              [Note: Used only for joining clinical tables + logging]
          (2)  ER_Status (see RECEPTOR_STATUSES list)
          (3)  PR_Status (see RECEPTOR_STATUSES list)
          (4)  HER2_Status (see RECEPTOR_STATUSES list)
          (5)  Histology_Type_HAOGUO (see HISTOLOGY_TYPES list)          [Note: When available, take precedent over CAMD's data]
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Check that the Biopsy_site_type of all samples from Hao's Dataset is 'Primary'
    assert (df["biopsy_site_type"] == 'Primary').all() == True

    required_columns = ['cBio_PatientID', 'cBio_SampleID', 'er', 'pr', 'her2___ihc',
                        'HER2_fish', 'histology']

    assert set(required_columns) <= set(df)
    df = df[required_columns]

    histology_colname = 'Histology_Type_HAOGUO'

    df.rename(columns={'cBio_PatientID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'cBio_SampleID': 'Tumor_Sample_Barcode'}, inplace=True)
    df.rename(columns={'er': 'ER_Status'}, inplace=True)
    df.rename(columns={'pr': 'PR_Status'}, inplace=True)
    df.rename(columns={'her2___ihc': 'HER2_IHC'}, inplace=True)
    df.rename(columns={'HER2_fish': 'HER2_FISH'}, inplace=True)
    df.rename(columns={'histology': histology_colname}, inplace=True)

    df = infer_HER2_status(df)
    df = df.drop(['HER2_IHC', 'HER2_FISH'], 1)

    # Fill-in missing values
    for colname in RECEPTOR_STATUS_COLS + [histology_colname]:
        df[colname] = df[colname].fillna(value='Unknown')

    # Standardize the values in the [ER,PR]_Status columns.
    for colname in ['ER_Status', 'PR_Status']:
        df[colname].replace('0.0', 'Negative', inplace=True)
        df[colname].replace('1.0', 'Positive', inplace=True)
        df[colname].replace('-1.0', 'Unknown', inplace=True)
        df[colname].replace('93.0', 'Unknown', inplace=True)

        # Map '2.0=Positive low (1-10)' to 'Negative' to be consistent with
        # Receptor Status data from public datasets (e.g. TCGA-BRCA).
        df[colname].replace('2.0', 'Negative', inplace=True)

    # Standardize the values in the Histology_Type column
    #    (1) Used the following source to guide the mapping:
    #           (A) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3047091/ (see Figure 1)
    #           (B) http://breast-cancer.ca/  (navigate to Section 5)
    df[histology_colname].replace('0.0',  'Other_Breast_Cancer',                inplace=True)  # 'DCIS'
    df[histology_colname].replace('1.0',  'Invasive_Ductal_Carcinoma',          inplace=True)  # 'Invasive ductal'
    df[histology_colname].replace('2.0',  'Invasive_Lobular_Carcinoma',         inplace=True)  # 'Invasive lobular'
    df[histology_colname].replace('3.0',  'Other_Invasive_Breast_Carcinoma',    inplace=True)  # 'Tubular'
    df[histology_colname].replace('4.0',  'Other_Invasive_Breast_Carcinoma',    inplace=True)  # 'Mucinous'
    df[histology_colname].replace('5.0',  'Other_Invasive_Breast_Carcinoma',    inplace=True)  # 'Micropapillary'
    df[histology_colname].replace('6.0',  'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)  # 'Mixed (IDC&ILS)'
    df[histology_colname].replace('99.0', 'Unknown',                            inplace=True)  # 'Other'
    df[histology_colname].replace('-1.0', 'Unknown',                            inplace=True)  # 'Unknown'

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df[histology_colname].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there are no duplicated 'Patient_ID' or 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()
    assert not df.duplicated(subset=['Patient_ID']).any()

    # Check that for all rows, Exome_ID has Sample_ID as its prefix.
    for index, row in df.iterrows():
        assert row['Tumor_Sample_Barcode'].startswith(row['Patient_ID'] + '_S')

    return df

def infer_HER2_status(df):
    """Infer the sample's HER2_Status from 'HER2_IHC' and 'HER2_FISH' data.
    """

    df['HER2_IHC'] = df['HER2_IHC'].fillna(value='Unknown')
    df['HER2_FISH'] = df['HER2_FISH'].fillna(value='Unknown')

    # Figure out the HER2_Status of the sample based on HER2 IHC and FISH results
    df['HER2_IHC'] = df['HER2_IHC'].fillna(value='Unknown')
    df['HER2_IHC'].replace('0.0',  'Negative',  inplace=True)
    df['HER2_IHC'].replace('2.0',  'Equivocal', inplace=True)
    df['HER2_IHC'].replace('3.0',  'Positive',  inplace=True)
    df['HER2_IHC'].replace('-1.0', 'Unknown',   inplace=True)
    df['HER2_IHC'].replace('93.0', 'Unknown',   inplace=True)

    df['HER2_FISH'] = df['HER2_FISH'].fillna(value='Unknown')
    df['HER2_FISH'].replace('0.0',  'Negative',  inplace=True)
    df['HER2_FISH'].replace('2.0',  'Equivocal', inplace=True)
    df['HER2_FISH'].replace('1.0',  'Positive',  inplace=True)
    df['HER2_FISH'].replace('-1.0', 'Unknown',   inplace=True)
    df['HER2_FISH'].replace('93.0', 'Unknown',   inplace=True)

    for colname in ['HER2_IHC', 'HER2_FISH']:
        assert set(df[colname].unique()) == set(['Positive', 'Negative', 'Equivocal', 'Unknown'])

    df['HER2_Status'] = "Unknown"
    df.loc[(df.HER2_IHC=='Negative'), 'HER2_Status'] = 'Negative'
    df.loc[(df.HER2_FISH=='Negative'), 'HER2_Status'] = 'Negative'
    df.loc[(df.HER2_IHC=='Positive'), 'HER2_Status'] = 'Positive'  # Overide the Negative
    df.loc[(df.HER2_FISH=='Positive'), 'HER2_Status'] = 'Positive'  # Overide the Negative

    return df

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_camd_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Clinical Data from DFCI CAMD lab")

    parser.add_argument("--in_camd_gender", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Patients' Gender Data from DFCI CAMD lab")

    parser.add_argument("--in_haoguo_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Clinical Data file obtain from Hao Guo (BCB).")

    parser.add_argument("--details_mode", action="store_true", default=False,
                        help="Output extra rows/columns. For behavior, see header documentation.")

    parser.add_argument("--out_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to output Clinical Data.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_camd_clinical:", repr(options.in_camd_clinical)
    print "##   in_camd_gender:", repr(options.in_camd_gender)
    print "##   in_haoguo_clinical:", repr(options.in_haoguo_clinical)
    print "##   details_mode:", repr(options.details_mode)
    print "##   out_clinical:", repr(options.out_clinical)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
