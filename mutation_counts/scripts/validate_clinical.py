#!/usr/bin/env python2

"""
SYNOPSIS
    Script to validate the Individual Cancer Center's Clinical Data files.

NOTES

    Individual Clinical Data file should be validated against this script before
    being uploaded/registered on IntelCCC's cluster.

    If Validation SUCCEEDED, script will exit with return code = 0.

    If Validation FAILED, script will exit with return code = 99.

    If an unexpected ERROR occur, the script should exit with other non-zero
    return code.

EXAMPLES

    # Validate the 'DFCI' Clinical Data file.
    #     Expected Exit Code: 0
    ./validate_clinical.py --in_clinical ../examples/inputs/private/clinical/dfci.clinical_data.r5.tsv

    # Validate the 'METABRIC' Clinical Data file (missing OPTIONAL 'Gender' column)
    #     Expected Exit Code: 0
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/metabric.clinical_data.tsv

    # Validate the 'SangerWGS' Clinical Data file..
    #     Expected Exit Code: 0
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv

    # Validate the 'TCGA-BRCA' Clinical Data file.
    #     Expected Exit Code: 0
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/tcga.clinical_data.tsv

    # Validate the 'MSK-IMPACT' Clinical Data file.
    #     Expected Exit Code: 0
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/genie_msk.clinical_data.tsv

EDGE_CASES

    # Edge Case #1: Exercise code for syntax errors and ability to detect errors/warnings.
    #     Gracefully Exit with Return Code: 99
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/metabric.clinical_data.tsv --test_code; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #2: Validate the OLD 'METABRIC' Clinical Data file (missing REQUIRED) 'Biopsy_Site_Type'
                    and 'Histology_Type' columns)
    #     Gracefully Exit with Return Code: 99
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/metabric.clinical_data.old.mc_v0.5.tsv

    # Edge Case #3: Input Clinical Data file doesn't exist.
    #     Gracefully Exit with Return Code: 99
    ./validate_clinical.py --in_clinical ../examples/inputs/public/clinical/missing.maf; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #4: Pass in a MAF file. Test if code can gracefully handle weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_clinical.py --in_clinical ../examples/inputs/public/metabric/metabric.maf; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #5: Pass in a BED file. Test if code can gracefully handle weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_clinical.py --in_clinical ../examples/inputs/public/metabric.targetedIntervals.r1.bed; \
    echo -e "##\n## RETURN_CODE: $?"

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (April-2017)
"""

import sys
import platform

# Confirm that this script is being ran with a CPython v2.7.X interpreter.
is_correct_interpreter = True
if platform.python_implementation() != "CPython":
    is_correct_interpreter = False

if sys.version_info[0:2] != (2, 7):
    is_correct_interpreter = False

if not is_correct_interpreter:
    sys.stderr.write("Error: Unsupported Python Interpreter:")
    sys.stderr.write(" (%s, %s).\n" % (platform.python_implementation(), sys.version_info[0:3]))
    sys.stderr.write("\n")
    sys.stderr.write("Please call 'validate_maf.py' using a CPython v2.7.X Interpreter.\n")
    sys.exit(1)

import os
import argparse
import copy
import numpy as np
import pandas as pd

VERSION = 'mc_v0.9_dev'

ERRORS = list()
WARNINGS = list()

VALIDATION_FAIL_RETURN_CODE = 99

def main(options):

    in_df = import_clinical(options.in_clinical)

    check_column_header(in_df)

    check_missing_data(in_df)

    check_sample_uniqueness(in_df)

    check_receptor_status_values(in_df)

    check_biopsy_site_type_values(in_df)

    check_histology_type_values(in_df)

    check_gender_values(in_df)

    check_center(in_df)

    summarize_and_exit()

def import_clinical(infile):
    """Import the Clinical Data file as a pandas data-frame.

    Notes
    -----
    (1) This infile is a tab-delimited Clinical Data file. For example, see:

            https://tinyurl.com/ztjodyf
    """

    try:
        in_df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header = 0)
    except Exception as E:
        sys.stderr.write("##\n")
        sys.stderr.write("## ERROR: Fail to import Clinical Data file: %s\n" % repr(E))
        sys.exit(VALIDATION_FAIL_RETURN_CODE)

    return in_df

def check_column_header(in_df):
    """Check for missing or extra Clinical Data columns in header."""

    if options.test_code:
        in_df = copy.deepcopy(in_df)
        in_df.drop('Center', axis=1, inplace=True)
        in_df.drop('PR_Status', axis=1, inplace=True)
        in_df['An_Extra_Column'] = 'FOO'

    expected_columns = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status',
                        'HER2_Status', 'Biopsy_Site_Type', 'Histology_Type', 'Gender']

    optional_columns = ['Gender']

    columns = list(in_df)

    extra_columns = list(set(columns) - set(expected_columns))
    missing_columns = list(set(expected_columns) - set(columns))

    # Handle optional columns
    missing_columns = list(set(missing_columns) - set(optional_columns))

    if len(extra_columns) > 0:
        err_msg = "Extra column(s) in Clinical header: %s" % extra_columns
        ERRORS.append(err_msg)

    if len(missing_columns) > 0:
        err_msg = "Missing expected column(s) in Clinical header: %s" % missing_columns
        ERRORS.append(err_msg)

def check_missing_data(in_df):
    """Check for missing data in agreed 'non_null' columns (which at this point
    is every columns)."""

    non_null_columns = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status',
                        'HER2_Status']

    if options.test_code:
        in_df = copy.deepcopy(in_df)
        TEST_ROW = 0
        in_df.loc[in_df.index[TEST_ROW], "Center"] = np.nan; TEST_ROW+=1

        for column in non_null_columns:
            in_df.loc[in_df.index[TEST_ROW], column] = np.nan
        TEST_ROW+=1

        for column in non_null_columns:
            in_df.loc[in_df.index[TEST_ROW], column] = np.nan; TEST_ROW+=1

        for column in list(in_df):
            in_df.loc[in_df.index[TEST_ROW], column] = np.nan
        TEST_ROW+=1

    err_msg_list = list()

    for column in non_null_columns:

        if column not in list(in_df):
            # err_msg already produced by check_column_header().
            continue

        null_counts =  in_df[column].isnull().values.sum()
        if null_counts != 0:
            err_msg_list.append([column, null_counts])

    if len(err_msg_list) != 0:
        err_msg = "Missing data in column(s): ["
        err_msg += ", ".join(["%s(rows=%d)" % (repr(x[0]), x[1]) for x in err_msg_list])
        err_msg += "]"
        ERRORS.append(err_msg)

def check_sample_uniqueness(in_df):
    """Check for duplicated samples across rows."""

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        for index in range(10):
            base_sample = in_df.iloc[TEST_ROW]['Tumor_Sample_Barcode']
            for num_duplicates in range(2+index):
                in_df.loc[in_df.index[TEST_ROW], 'Tumor_Sample_Barcode'] = base_sample
                TEST_ROW += 1;

    sample_count = dict()

    if 'Tumor_Sample_Barcode' not in list(in_df):
        # err_msg already produced by check_column_header().
        return

    for index, row in in_df.iterrows():
        sample = row['Tumor_Sample_Barcode']

        if sample not in sample_count:
            sample_count[sample] = 0

        sample_count[sample] += 1

    duplicated_sample_strs = list()

    for sample, count in sample_count.iteritems():
        if count > 1:
            duplicated_sample_strs.append('%s (count=%s)' % (sample, count))

    duplicated_sample_strs = sorted(duplicated_sample_strs)

    if len(duplicated_sample_strs) > 0:
        max_show = 5
        pural = 's' if len(duplicated_sample_strs) > 1 else ''
        err_msg = "Found %s duplicated tumor_sample_barcode%s. E.g. %s" % (len(duplicated_sample_strs),
                                                                           pural,
                                                                           duplicated_sample_strs[:max_show])
        ERRORS.append(err_msg)

def check_receptor_status_values(in_df):
    """Check for invalid values in the receptor status columns."""

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        for colname in ['ER_Status', 'PR_Status', 'HER2_Status']:
            in_df.loc[in_df.index[TEST_ROW], colname] = 'no_data_supplied'; TEST_ROW+=1
            in_df.loc[in_df.index[TEST_ROW], colname] = 'Indeterminate'; TEST_ROW+=1
            in_df.loc[in_df.index[TEST_ROW], colname] = 'NA'; TEST_ROW+=1
            in_df.loc[in_df.index[TEST_ROW], colname] = '.'; TEST_ROW+=1
            in_df.loc[in_df.index[TEST_ROW], colname] = ''; TEST_ROW+=1

    for colname in ['ER_Status', 'PR_Status', 'HER2_Status']:
        valid_values = ['Positive', 'Negative', 'Unknown']

        if colname not in list(in_df):
            # err_msg already produced by check_column_header().
            return

        observed_values = list(in_df[colname].unique())

        invalid_values = set(observed_values) - set(valid_values)

        if len(invalid_values) > 0:
            err_msg = "Invalid value(s) in '%s' column: %s" % (colname, list(invalid_values))
            ERRORS.append(err_msg)

def check_biopsy_site_type_values(in_df):
    """Check for invalid values in the 'Biopsy_Site_Type' column"""

    colname = 'Biopsy_Site_Type'

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Local Recurrence'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Metastatic Recurrence'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Metastasis'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Any/Other'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'no_data_supplied'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Unspecified'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Not Applicable'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'NA'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = '.'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = ''; TEST_ROW+=1

    if colname not in list(in_df):
        # 'Biopsy_Site_Type' is an optional column.
        return

    valid_values = ['Primary', 'Local_Recurrence', 'Metastatic', 'Unknown']

    observed_values = list(in_df[colname].unique())

    invalid_values = set(observed_values) - set(valid_values)

    if len(invalid_values) > 0:
        err_msg = "Invalid value(s) in '%s' column: %s" % (colname, list(invalid_values))
        ERRORS.append(err_msg)

def check_histology_type_values(in_df):
    """Check for invalid values in the 'Histology_Type' column"""

    colname = 'Histology_Type'

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        in_df.loc[in_df.index[TEST_ROW], colname] = 'invasive ductal'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'invasive lobular'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'invasive mixed'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'IDC'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'ILC'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'no_data_supplied'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Unspecified'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Not Applicable'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'NA'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = '.'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = ''; TEST_ROW+=1

    if colname not in list(in_df):
        # 'Biopsy_Site_Type' is an optional column.
        return

    valid_values = ['Invasive_Ductal_Carcinoma', 'Invasive_Lobular_Carcinoma', 'Mixed_Ductal_and_Lobular_Carcinoma',
                    'Other_Invasive_Breast_Carcinoma', 'Other_Breast_Cancer', 'Unknown']

    observed_values = list(in_df[colname].unique())

    invalid_values = set(observed_values) - set(valid_values)

    if len(invalid_values) > 0:
        err_msg = "Invalid value(s) in '%s' column: %s" % (colname, list(invalid_values))
        ERRORS.append(err_msg)

def check_gender_values(in_df):
    """Check for invalid values in the 'Gender' column"""

    colname = 'Gender'

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        in_df.loc[in_df.index[TEST_ROW], colname] = 'F'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'M'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'female'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'male'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'no_data_supplied'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Unspecified'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'Not Applicable'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = 'NA'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = '.'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], colname] = ''; TEST_ROW+=1

    if colname not in list(in_df):
        # 'Gender' is an optional column.
        return

    valid_values = ['Male', 'Female', 'Unknown']

    observed_values = list(in_df[colname].unique())

    invalid_values = set(observed_values) - set(valid_values)

    if len(invalid_values) > 0:
        err_msg = "Invalid value(s) in '%s' column: %s" % (colname, list(invalid_values))
        ERRORS.append(err_msg)

def check_center(in_df):
    """Check for consistency in the 'Center' value in the rows."""

    if options.test_code:
        TEST_ROW = 0
        in_df = copy.deepcopy(in_df)
        in_df.loc[in_df.index[TEST_ROW], 'Center'] = 'Foo'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], 'Center'] = 'Bar'; TEST_ROW+=1
        in_df.loc[in_df.index[TEST_ROW], 'Center'] = ''; TEST_ROW+=1

    if 'Center' not in list(in_df):
        # err_msg already produced by check_column_header().
        return

    values = list(in_df['Center'].unique())

    if len(values) > 1:
        err_msg = "Expect single value in 'Center' column across all rows. Found multiple values in 'Center' column: %s" % list(values)
        ERRORS.append(err_msg)

def summarize_and_exit():
    """Summarize the MAF validation result and exit (with status 1 if found
    errors)."""

    if len(ERRORS) != 0:
        sys.stderr.write("##\n")
        sys.stderr.write("## Clinical Data Validation FAILED. Please address the following critical ERROR(s):\n")
        for index, msg in enumerate(ERRORS, start=1):
            sys.stderr.write("##\n")
            sys.stderr.write("##    %d. %s\n" % (index, msg))


        if len(WARNINGS) != 0:
            sys.stderr.write("##\n")
            sys.stderr.write("## Additionally encountered the following WARNING(s):\n")
            for index, msg in enumerate(WARNINGS, start=1):
                sys.stderr.write("##\n")
                sys.stderr.write("##    %d. %s\n" % (index, msg))

        sys.exit(VALIDATION_FAIL_RETURN_CODE)

    elif len(WARNINGS) != 0:
        sys.stderr.write("##\n")
        sys.stderr.write("## Clinical Data Validation SUCCEEDED with the following WARNING(s):\n")
        for index, msg in enumerate(WARNINGS, start=1):
            sys.stderr.write("##\n")
            sys.stderr.write("##    %d. %s\n" % (index, msg))

        sys.exit(0)

    else:
        sys.stderr.write("##\n")
        sys.stderr.write("## Clinical Data Validation SUCCEEDED with no WARNING.\n")
        sys.exit(0)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input Clinical data file for validation.")

    parser.add_argument("--test_code", action='store_true',
                        help="Exercise code for syntax errors and ability to detect errors/warnings.")

    options = parser.parse_args()

    print "## Program Version:", repr(VERSION)
    print "## Validating Clinical Data file %s." % repr(options.in_clinical)

    main(options)

    sys.exit(0)
