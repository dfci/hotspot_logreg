#!/usr/bin/env python2

"""
SYNOPSIS
    Merge data from multiple input Breast Cancer Clinical files into a single
    output Clinical data file.

NOTES

    (a) Check that all input clinical files have the expected columns.

    (b) Append 'Center' name to 'Tumor_Sample_Barcode' as a way to help prevent
        samples from different input MAF files from colliding. Motivation is to
        ensure that 'Tumor_Sample_Barcode' can be used on a unique identifier of
        samples in downstream analysis steps.

        - After modifing the 'Tumor_Sample_Barcode' values, then check that
          there is no collision between the 'Tumor_Sample_Barcode' values of
          the different input MAF files.

EXAMPLES

    # Merge clinical data from DFCI, MSK-IMPACT, METABRIC, SangerWGS, and TCGA
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/private/clinical/dfci.clinical_data.r5.tsv \
                       ../examples/inputs/public/clinical/genie_msk.clinical_data.tsv \
                       ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
                       ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv \
                       ../examples/inputs/public/clinical/tcga.clinical_data.tsv \
        --out_clinical private/merged.clinical_data.r13.tsv

    # Merge clinical data from only Public sources (MSK-IMPACT, METABRIC, SangerWGS, and TCGA)
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/genie_msk.clinical_data.tsv \
                       ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
                       ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv \
                       ../examples/inputs/public/clinical/tcga.clinical_data.tsv \
        --out_clinical private/merged.clinical_data.only_public.r13.tsv

    # Merge clinical data from DFCI, MSK-IMPACT, METABRIC, SangerWGS, and TCGA.
    # However for the METABRIC samples, include only the ER+ subset.
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/private/clinical/dfci.clinical_data.r5.tsv \
                       ../examples/inputs/public/clinical/genie_msk.clinical_data.tsv \
                       ../examples/inputs/public/clinical/metabric.clinical_data.er_pos_subset.tsv \
                       ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv \
                       ../examples/inputs/public/clinical/tcga.clinical_data.tsv \
        --out_clinical private/merged.clinical_data.er_pos_subset_metabric.r13.tsv

EDGE_CASES

    Edge Case #1: Single input MAF file.
        Expected Behavior: Gracefully handle and generate correct output.
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
        --out_clinical private/tmp.clinical_data.tsv

    Edge Case #2: One of the clinical dataset (METABRIC) is missing 'Gender' column
        Expected Behavior: Gracefully handle and fill in default values for the missing columns
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
                       ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv \
        --out_clinical private/tmp.clinical_data.tsv

    Edge Case #3: One of the clinical dataset (METABRIC) is missing 'Biopsy_Site_Type'
                  and 'Histology_Type' columns.
        Expected Behavior: Raise an error, since these column are no longer 'optional'.
    ./merge_clinical.py \
         --in_clinicals ../examples/inputs/public/clinical/metabric.clinical_data.old.mc_v0.5.tsv \
                        ../examples/inputs/public/clinical/sanger_wgs.clinical_data.tsv \
         --out_clinical private/tmp.clinical_data.tsv

    Edge Case #4: Some of the input MAF files does not exist.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
                       ../examples/inputs/public/clinical/missing.maf \
        --out_clinical private/tmp.clinical_data.tsv

    Edge Case #5: Collision between basenames of the different input_clinical files.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
                       ../examples/inputs/public/clinical/metabric.clinical_data.tsv \
        --out_clinical private/tmp.clinical_data.tsv

    Edge Case #6: Collision between the 'Tumor_Sample_Barcode' values of the different input clinical files.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_clinical.py \
        --in_clinicals ../examples/inputs/public/clinical/collison_1.clinical_data.tsv \
                       ../examples/inputs/public/clinical/collison_2.clinical_data.tsv \
        --out_clinical private/tmp.clinical_data.tsv

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
    sys.stderr.write("Please call 'merge_maf.py' using a CPython v2.7.X Interpreter.\n")
    sys.exit(1)

import os
import time
import argparse
import pandas as pd

VERSION = 'mc_v0.9_dev'

expected_columns = ['Tumor_Sample_Barcode', 'Center', 'ER_Status', 'PR_Status', 'HER2_Status', 'Biopsy_Site_Type', 'Histology_Type', 'Gender']
optional_columns = ['Gender']

def main(options):

    in_clinical_dfs = list()
    in_centers = set()  # For 'Center' collision check
    in_center2samples = dict()  # For 'Sample' collision check

    for in_clinical_file in options.in_clinicals:

        in_df = pd.read_table(in_clinical_file, sep="\t", dtype=str, comment="#", header=0)

        # Ensure that there is no missing data in any of the columns.
        assert not in_df.isnull().values.any()

        # Ensure that all the samples in each input clinical data file have the
        # same center name.
        assert in_df['Center'].nunique() == 1

        center = list(in_df['Center'].unique())[0]

        # Ensure that there all no duplicate centers.
        if center in in_centers:
            raise Exception("Duplicated Center '%s' across input clinical files." % center)
        in_centers.add(center)

        # If an optional column is missing, then create and fill it will the
        # default value.
        for column in optional_columns:
            if column not in set(in_df):
                if column == 'Gender':
                    print "##"
                    print "## WARNING: Fill-in missing 'Gender' column for center '%s'." % center
                    in_df[column] = 'Unknown'
                else:
                    raise Exception("Invalid optional column '%s'" % column)

        columns = set(in_df)
        if columns != set(expected_columns):
            extra_columns = list(set(columns) - set(expected_columns))
            missing_columns = list(set(expected_columns) - set(columns))

            # Handle optional columns
            missing_columns = list(set(missing_columns) - set(optional_columns))

            if len(extra_columns) > 0:
                raise Exception("Found extra clinical column(s): %s" % extra_columns)

            if len(missing_columns) > 0:
                raise Exception("Missing expected clinical column(s): %s" % missing_columns)

        # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values with center.
        assert not in_df.duplicated(subset=['Tumor_Sample_Barcode']).any()

        # Append 'Center' name to 'Tumor_Sample_Barcode' as a way to help prevent
        # sample_names from different centers from colliding.
        in_df["Tumor_Sample_Barcode"] = center + "_" + in_df["Tumor_Sample_Barcode"]

        check_tumor_sample_barcode_collisions(in_df, center, in_center2samples)

        # Keep only the desired columns in the desired order.
        in_df = in_df[expected_columns]

        in_clinical_dfs.append(in_df)

    # Concatenate and output the Clinical Data
    out_df = pd.concat(in_clinical_dfs, axis=0)

    # Ensure that there is no missing data in any of the columns.
    assert not out_df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values across
    # centers (note: this is the modified 'Tumor_Sample_Barcode').
    assert not out_df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    out_df.to_csv(options.out_clinical, sep="\t", na_rep='', index=False)

    # Print Merged Clinical Data Summary for debugging/logging purposes.
    print "##", "-" * 50
    print "## Merged Clinical Data Summary:"
    print "##   Total Num Samples:", format(len(out_df['Tumor_Sample_Barcode'].unique()), ',d')
    print "##", "-" * 50

def check_tumor_sample_barcode_collisions(curr_df, curr_center, prev_center2samples):
    """Ensure no collision between 'Tumor_Sample_Barcode' values from 'curr_df'
    and the previously processed clinical files ('prev_center2samples').

    Notes
    -----
    (1) If there is collision, this function will raise an Exception.

    (2) If there is no collision, this function will add the samples from
        'curr_maf_df' to 'prev_maf2samples' dict and then return.
    """

    samples = set(curr_df["Tumor_Sample_Barcode"].unique())
    assert curr_center not in prev_center2samples

    for prev_center, prev_samples in prev_center2samples.iteritems():
        duplicates = sorted(list(prev_samples.intersection(samples)))

        if len(duplicates) != 0:
            max_show = max(len(duplicates), 5)
            err_msg = "Found %s duplicated 'Tumor_Sample_Barcode' values between" % len(duplicates)
            err_msg += " center_1 %s and center_2 %s." % (repr(curr_center), repr(prev_center))
            if len(duplicates) > 5:
                err_msg += " Here is a subset of the duplicated values: %s" % duplicates[:5]
            else:
                err_msg += " Here are the duplicated values: %s" % duplicates

            raise Exception(err_msg)

    prev_center2samples[curr_center] = samples

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_clinicals", action="store", required=True, nargs='+',
                        metavar='FILE',
                        help="List of paths to input Clinical Data files.")

    parser.add_argument("--out_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output (merged) Clinical Data file.")

    options = parser.parse_args()

    for index in range(len(options.in_clinicals)):
        options.in_clinicals[index] = os.path.abspath(options.in_clinicals[index])

    options.out_clinical= os.path.abspath(options.out_clinical)

    print "##", "-" * 50
    print "## Program Version:", repr(VERSION)
    print "## Specified Options:"
    print "##   in_clinicals:"
    for index, in_clinical in enumerate(options.in_clinicals, start=1):
        print "##       (%d) %s" % (index, repr(in_clinical))
    print "##   out_clinical:", repr(options.out_clinical)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
