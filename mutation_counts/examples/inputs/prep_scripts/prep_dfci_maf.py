#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the DFCI Profile Breast Cancer data into a MAF format that is suitable
    for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Import the 'details' preprocessed clinical data file. This clinical data
        contains:
            (A) A 'Panel_Version' column which indicates the panel_version of each
                sample.

            (B) A 'Exclude_Sample' column which indicate whether the sample
                would be excluded (filtered-out) from the outfile based on
                clinical filtering citeria.

                Note
                ----
                See 'prep_dfci_clinical.py' script for more information regarding
                the clinical filtering citeria to exclude samples.

    (2) Assume that REF/ALT/COORDINATES values in input ('--in_genomics') already
        satisfy the MAF-specification.

    (3) There are currently three ONCOPANEL versions ('1', '2', '3')

EXAMPLES

    # Extract genomics data into suitable output MAF format. One outfile per
    # panel (v1, v2, v3)
    ./prep_dfci_maf.py \
        --in_genomics ../private/raw/06Feb2017/vep_annotated.maf \
        --in_clinical ../private/clinical/dfci.clinical_data.r5.details.tsv \
        --panel_versions 1 2 3 \
        --out_mafs ../private/dfci.DFCI-ONCOPANEL-1.r4.check.maf \
                   ../private/dfci.DFCI-ONCOPANEL-2.r4.check.maf \
                   ../private/dfci.DFCI-ONCOPANEL-3.r4.check.maf

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Mar-2017)
"""


import sys
import os
import time
import copy
import subprocess
import argparse

import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

MAF_columns = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
               "Chromosome", "Start_Position", "End_Position", "Strand",
               "Variant_Classification", "Variant_Type", "Reference_Allele",
               "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS",
               "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
               "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
               "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
               "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
               "Verification_Status", "Validation_Status", "Mutation_Status",
               "Sequencing_Phase", "Sequence_Source", "Validation_Method",
               "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID",
               "Matched_Norm_Sample_UUID"]

PANEL_VERSIONS = ['1', '2', '3']

def main(options):

    genomics_df = import_dfci_genomics(options.in_genomics)
    clinical_df = import_prepped_clinical(options.in_clinical)

    # Check that each sample in genomics_df has corresponding sample in clinical_df.
    # However, note that the converse is not true.
    assert set(genomics_df['Tumor_Sample_Barcode'].unique()) <= set(clinical_df['Tumor_Sample_Barcode'].unique())

    # Filter for data with the desired Panel(s)
    clinical_df = clinical_df[clinical_df['Panel_Version'].isin(options.panel_versions)]

    for index in range(len(options.panel_versions)):

        panel_version = options.panel_versions[index]
        out_maf = options.out_mafs[index]

        # Filter for match Panel_Version
        filtered_clinical_df = clinical_df[clinical_df['Panel_Version'] == panel_version]

        # Filter-out 'Exclude_Sample' (see 'prep_dfci_clinical.py' for details).
        filtered_clinical_df = filtered_clinical_df[filtered_clinical_df['Exclude_Sample'] == 'False']

        # Merge the filtered clinical and genomics df.
        out_df = pd.merge(genomics_df,
                          filtered_clinical_df,
                          how="inner",  # INNER JOIN
                          left_on="Tumor_Sample_Barcode",
                          right_on="Tumor_Sample_Barcode",
                          sort=False,
                          indicator='indicator_column')

        # Consistency checks
        # (1) This should be true, since using INNER JOIN.
        assert set(out_df['indicator_column'].unique()) == set(['both'])

        # (2) Check for duplicate rows
        any_duplicates = out_df.duplicated(subset=["Tumor_Sample_Barcode", "Start_Position",
                                                   "End_Position", "Reference_Allele",
                                                   "Tumor_Seq_Allele2"]).any()

        assert any_duplicates == False

        # (3) Check that the values in these two columns are the same for every row.
        assert (out_df['Tumor_Seq_Allele1'] == out_df['Reference_Allele']).all()

        # Set Center to 'DFCI'
        out_df['Center'] = 'DFCI'

        # If a standard MAF column doesn't exist, then add column with null ('') values.
        existing_columns = set(out_df.columns.values)
        for colname in MAF_columns:
            if colname not in existing_columns:
                out_df[colname] = ''

        # Assert that there are valid data-rows in the out_df
        assert len(out_df) > 0

        # Reorder and keep only desired columns
        out_df = out_df[MAF_columns]

        fout = open(out_maf, 'w')
        fout.write("#version 2.4\n")
        out_df.to_csv(fout, sep="\t", na_rep='', index=False)
        fout.close()

        print "##", "-" * 50
        print "## Outfile Summary (Panel: 'OncoPanel-v%s'):" % panel_version
        print "##   Total # Samples (Clinical):", filtered_clinical_df['Tumor_Sample_Barcode'].nunique()
        print "##   Total # Samples (Genomics):", out_df['Tumor_Sample_Barcode'].nunique()
        print "##   Total # Variant Calls:", len(out_df)
        print "##", "-" * 50

def import_prepped_clinical(infile):
    """Import the prepped clinical file.

    Notes
    -----
    This infile should be a clinical data file that has been preprocessed by
    the correspond 'prep_*_clinical.py' script using the '--details_mode'.

    Extract the following columns from the input dfci clinical file:
          (1)  Tumor_Sample_Barcode (e.g. 'CBIO_P10001_S1')
          (2)  Panel_Version (see PANEL_VERSIONS list)
          (3)  Exclude_Sample ('True', 'False')
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Check that all the require columns exist
    required_columns = ['Tumor_Sample_Barcode', 'Panel_Version', 'Exclude_Sample']

    assert set(required_columns) <= set(df)
    df = df[required_columns]

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Ensure that there is no unexpected 'Panel_Version' values.
    assert set(df['Panel_Version'].unique()) <= set(PANEL_VERSIONS)

    # Ensure that there is no unexpected 'Exclude_Sample' values.
    assert set(df['Exclude_Sample'].unique()) <= set(['True', 'False'])

    return df

def import_dfci_genomics(infile):
    """Import the DFCI (cBioOne) genomics file (which actually is a MAF file
    itself).

    Notes
    -----
    Extract the following columns from the input dfci genomics file:
        (1) Tumor_Sample_Barcode
        (2) Chromosome
        (3) Start_Position
        (4) End_Position
        (5) Strand
        (6) Reference_Allele
        (7) Tumor_Seq_Allele1
        (8) Tumor_Seq_Allele2
        (9) NCBI_Build
    """

    df = pd.read_table(infile, sep="\t", dtype=str, header=0)

    keep_colnames = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                     "End_Position", "Strand", "Reference_Allele",
                     "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "NCBI_Build"]

    # Keep only the columns we need.
    df = df[keep_colnames]

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there is no unexpected value of Strand column.
    assert df["Strand"].unique() == ["+"]

    # Ensure that there is no unexpected value of NCBI_Build column.
    assert df["NCBI_Build"].unique() == ["GRCh37"]

    return df

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_genomics", action="store", required=True,
                        metavar='FILE',
                        help="Path to input DFCI Genomics Data.")

    parser.add_argument("--in_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to input (Detailed Proprocessed) DFCI Clinical Data.")

    parser.add_argument("--panel_versions", action="store", required=True,
                         nargs='+',
                        help="OncoPanel Version filter.")

    parser.add_argument("--out_mafs", action="store", required=True,
                         nargs='+',
                         help="Path to the output MAF file(s).")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_genomics: ", repr(options.in_genomics)
    print "##   in_clinical: ", repr(options.in_clinical)
    print "##   panel_versions: ", repr(options.panel_versions)
    print "##   out_mafs:", repr(options.out_mafs)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
