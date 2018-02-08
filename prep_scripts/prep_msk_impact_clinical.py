#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the MSK-IMPACT Breast Cancer Clinical data into a format that is
    suitable for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Support 'detail' mode, which will:

          (A) Include 'Panel_Version' column (note: this is mapped from 'SEQ_ASSAY_ID')

          (B) Retain samples which will be normally filtered out (e.g. duplicated
              samples from a patient, Breast Sarcoma samples) and instead
              add 'Exclude_Sample' column  ('True', 'False') and 'Exclude_Reason'
              column (text describing the reason to exclude sample).

              For sample exclusion logic, see filter_clinical_df() function.

        Note
        ----
        Purpose of the 'detail' mode is to include extra rows/columns in the
        preprocessed clinical outfile that can be consumed by the 'prep_msk_impact_maf.py'
        script BUT is not used/supported by the actual mutation  counts pipeline.

    (2) According to GENIE data guide (http://www.aacr.org/Documents/GENIEDataGuide.pdf):

          - Clinical Data's 'CANCER_TYPE' column corresponds to:
              The primary cancer diagnosis "main type", based on the
              OncoTree ontology (http://cbioportal.org/oncotree).

          - Clinical Data's 'CANCER_TYPE_DETAILED' column corresponds to:
              The primary cancer diagnosis label, based on the OncoTree ontology
              (http://cbioportal.org/oncotree).

    (3) There seem be two different MSK-IMPACT Panels:

          - 'MSK-IMPACT341' (Version 1 of the panel, containing 341 Genes)
          - 'MSK-IMPACT410' (Version 2 of the panel, containing 410 Genes)

EXAMPLES

    # (1) Generate the 'regular' preprocessed MSK-IMPACT clinical outfile.
    ./prep_msk_impact_clinical.py \
        --in_clinical ../public/genie/raw/03Feb2017/data_clinical.txt \
        --out_clinical ../public/clinical/genie_msk.clinical_data.tsv

    # (2) Generate the 'detail' preprocessed MSK-IMPACT clinical outfile.
    ./prep_msk_impact_clinical.py \
        --in_clinical ../public/genie/raw/03Feb2017/data_clinical.txt \
        --details_mode \
        --out_clinical ../public/clinical/genie_msk.clinical_data.details.tsv


AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Mar-2017)
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

BIOPSY_SITE_TYPES = ['Primary', 'Local_Recurrence', 'Metastatic', 'Unknown']

HISTOLOGY_TYPES = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                   'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

GENDER_TYPES = ['Female', 'Male', 'Unknown']

PANEL_VERSIONS = ["MSK-IMPACT341", "MSK-IMPACT410"]

def main(options):

    clinical_df = import_genie_clinical(options.in_clinical)

    # For debugging/logging purposes
    print_patients_with_inconsistent_cancer_types(clinical_df)

    clinical_df = filter_clinical_df(clinical_df, options.details_mode)

    # Exclude 'Breast Sarcoma' (only 2 counts in GENIE MSK-IMPACT Breast Cancer dataset)
    # clinical_df = clinical_df[clinical_df['CANCER_TYPE'] == 'Breast Cancer']

    # Set Receptor status columns to 'Unknown' since this information is not
    # available.
    clinical_df['ER_Status'] = 'Unknown'
    clinical_df['PR_Status'] = 'Unknown'
    clinical_df['HER2_Status'] = 'Unknown'

    # Set Center to 'GENIE-MSK'
    clinical_df['Center'] = 'GENIE-MSK'

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

         (2) Exclude Breast_Sarcoma samples (2 in MSK dataset)
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

        age_str = row['Age_At_Seq_Report']
        if age_str == "<18y":
            age = 17.0
        elif age_str == ">89y":
            age = 90.0
        else:
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

def import_genie_clinical(infile):
    """Import the MSK-IMPACT clinical file (obtained from GENIE project)

    Notes
    -----
    - Extract the following data:
          (1)  Patient_ID (e.g. 'CBIO_P10001')                   [Note: Used only debugging/logging]
          (2)  Tumor_Sample_Barcode (e.g. 'CBIO_P10001_S1')
          (3)  Biopsy_Site_Type (see BIOPSY_SITE_TYPES list)
          (4)  Histology_Type (see HISTOLOGY_TYPES list)
          (5)  Gender (see GENDER_TYPES list)
          (6)  Panel_Version (see PANEL_VERSIONS list)           [Note: Only output in detailed mode]
          (7)  CANCER_TYPE                                       [Note: Only for debugging/logging purposes]
          (8)  Is_Breast_Sarcoma (True, False)                   [Note: Only for downstream filtering purposes]
          (9)  Is_Breast_Carcinoma (True, False)                 [Note: Only for downstream filtering purposes]
          (0)  Age_At_Seq_Report                                 [Note: Only for downstream filtering purposes]
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Filter for only data from 'MSK' Center.
    df = df[df["CENTER"] == "MSK"]

    df['Is_Breast_Sarcoma'] = (df['CANCER_TYPE'] == 'Breast Sarcoma')  # Not consider a 'Breast Cancer'
    df['Is_Breast_Carcinoma'] = (df['CANCER_TYPE'] == 'Breast Cancer')

    # Retain only the required columns
    required_columns = ['PATIENT_ID', 'SAMPLE_ID', 'SAMPLE_TYPE', 'CANCER_TYPE_DETAILED', 'SEX', 'SEQ_ASSAY_ID', 'CANCER_TYPE',
                        'Is_Breast_Carcinoma', 'Is_Breast_Sarcoma', 'AGE_AT_SEQ_REPORT']

    assert set(required_columns) <= set(df)
    df = df[required_columns]

    df.rename(columns={'PATIENT_ID': 'Patient_ID'}, inplace=True)
    df.rename(columns={'SAMPLE_ID': 'Tumor_Sample_Barcode'}, inplace=True)
    df.rename(columns={'SAMPLE_TYPE': 'Biopsy_Site_Type'}, inplace=True)
    df.rename(columns={'CANCER_TYPE_DETAILED': 'Histology_Type'}, inplace=True)
    df.rename(columns={'SEX': 'Gender'}, inplace=True)
    df.rename(columns={'SEQ_ASSAY_ID': 'Panel_Version'}, inplace=True)
    df.rename(columns={'AGE_AT_SEQ_REPORT': 'Age_At_Seq_Report'}, inplace=True)

    # Standardize the values in the 'Biopsy_Site_type' column
    df['Biopsy_Site_Type'].replace('Metastasis', 'Metastatic', inplace=True)

    # Standardize the values in the 'Histology_Type' column
    df['Histology_Type'].replace('Breast Invasive Ductal Carcinoma',          'Invasive_Ductal_Carcinoma',          inplace=True)
    df['Histology_Type'].replace('Breast Invasive Lobular Carcinoma',         'Invasive_Lobular_Carcinoma',         inplace=True)
    df['Histology_Type'].replace('Breast Mixed Ductal and Lobular Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma', inplace=True)
    df['Histology_Type'].replace('Invasive Breast Carcinoma',                 'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Adenoid Cystic Breast Cancer',              'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Breast Invasive Carcinoma, NOS',            'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Breast Invasive Cancer, NOS',               'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Breast Invasive Mixed Mucinous Carcinoma',  'Other_Invasive_Breast_Carcinoma',    inplace=True)
    df['Histology_Type'].replace('Metaplastic Breast Cancer',                 'Other_Breast_Cancer',                inplace=True)
    df['Histology_Type'].replace('Metaplastic Carcinosarcoma',                'Other_Breast_Cancer',                inplace=True)
    df['Histology_Type'].replace('Malignant Phyllodes Tumor of the Breast',   'Other_Breast_Cancer',                inplace=True)

    df.loc[~df.Is_Breast_Carcinoma, 'Histology_Type'] = 'Not_Breast_Carcinoma'

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Check that for all rows, Exome_ID has Sample_ID as its prefix.
    for index, row in df.iterrows():
        assert row['Tumor_Sample_Barcode'].startswith(row['Patient_ID'] + '-T')

    # Ensure that there is no unexpected 'Biopsy_Site_Type' values.
    assert set(df['Biopsy_Site_Type'].unique()) <= set(BIOPSY_SITE_TYPES)

    # Ensure that there is no unexpected 'Histology_Type' values.
    assert set(df[df['Is_Breast_Carcinoma']]['Histology_Type'].unique()) <= set(HISTOLOGY_TYPES)

    # Ensure that there is no unexpected 'Gender' values.
    assert set(df["Gender"].unique()) <= set(GENDER_TYPES)

    # Ensure that there is no unexpected 'Panel_Version' values.
    assert set(df["Panel_Version"].unique()) <= set(PANEL_VERSIONS)

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
