#!/usr/bin/env python2

"""
SYNOPSIS
    Prepare the TCGA-BRCA data into a MAF format that is suitable for input into
    the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Assume that "Reference_Allele"/"Tumor_Seq_Allele1"/"Tumor_Seq_Allele2"
        values in input ('--in_genomics') already satisfy the MAF-specification.

EXAMPLES

    # Extract genomics data into suitable output MAF format.
    ./prep_tcga_maf.py \
        --in_genomics ../public/tcga_brca/raw/brca_tcga_pub2015/brca_tcga_pub2015/data_mutations_extended.txt \
        --out_maf ../public/tcga_brca/tcga.cell_2015.r2.check.maf

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (April-2017)
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

def main(options):

    out_df = import_cbioportal_genomics(options.in_genomics)

    out_df = update_tumor_seq_allele1(out_df)
    out_df = filter_non_primary_chromosomes(out_df)

    # Check for duplicate rows
    any_duplicates = out_df.duplicated(subset=["Tumor_Sample_Barcode", "Start_Position",
                                               "End_Position", "Reference_Allele",
                                               "Tumor_Seq_Allele2"]).any()

    assert any_duplicates == False

    # Check that the values in these two columns are the same for every row.
    assert (out_df['Tumor_Seq_Allele1'] == out_df['Reference_Allele']).all()

    # Set Center to 'TCGA'
    out_df['Center'] = 'TCGA'

    # If a standard MAF column doesn't exist, then add column with null ('') values.
    existing_columns = set(out_df.columns.values)
    for colname in MAF_columns:
        if colname not in existing_columns:
            out_df[colname] = ''

    # Reorder and keep only desired columns
    out_df = out_df[MAF_columns]

    fout = open(options.out_maf, 'w')
    fout.write("#version 2.4\n")
    out_df.to_csv(fout, sep="\t", na_rep='', index=False)
    fout.close()

    print "##", "-" * 50
    print "## Outfile Summary:"
    print "##   Total # Samples (Genomics):", out_df['Tumor_Sample_Barcode'].nunique()
    print "##   Total # Variant Calls:", len(out_df)
    print "##", "-" * 50

def import_cbioportal_genomics(infile):
    """Import the TCGA-BRCA genomics file (which actually is itself a MAF file
    downloaded from cBioPortal).

    Notes
    -----
    Extract the following columns from the input genomics file:
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

    df = pd.read_table(infile, sep="\t", dtype=str, header=0, comment="#")

    keep_colnames = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                     "End_Position", "Strand", "Reference_Allele",
                     "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "NCBI_Build"]

    # Keep only the columns we need.
    df = df[keep_colnames]

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    invalid_values = list(set(df['Strand'].unique()) - set(['+']))
    if len(invalid_values) > 0:
        print ""
        print "## WARNING: Found invalid 'Strand' values:", invalid_values
        print "## WARNING: Set all 'Strand' values to '+' for now, but may want to revisit/investigate this issue."
        df['Strand'] = '+'

    invalid_values = list(set(df['NCBI_Build'].unique()) - set(['GRCh37']))
    if len(invalid_values) > 0:
        print "##"
        print "## WARNING: Found invalid 'NCBI_Build' values:", invalid_values
        print "## WARNING: Set all 'NCBI_Build' values to 'GRCh37' for now."
        df['NCBI_Build'] = 'GRCh37'

    return df

def update_tumor_seq_allele1(df):
    """Check and update the value in the 'Tumor_Seq_Allele1' in the cases where
    'Tumor_Seq_Allele1' does not equal 'Reference_Allele'

    Notes
    -----
    (1) Make sure that in these cases, 'Tumor_Seq_Allele1' equals 'Tumor_Seq_Allele2'

    (2) Then reset 'Tumor_Seq_Allele1' to equal 'Reference_Allele'
    """

    indices = df.index[df['Tumor_Seq_Allele1'] != df['Reference_Allele']]

    for index in indices:

        reference_allele = df.get_value(index, 'Reference_Allele')
        tumor_seq_allele1 = df.get_value(index, 'Tumor_Seq_Allele1')
        tumor_seq_allele2 = df.get_value(index, 'Tumor_Seq_Allele2')

        assert tumor_seq_allele1 == tumor_seq_allele2

        df.set_value(index, 'Tumor_Seq_Allele1', reference_allele)

    return df

def filter_non_primary_chromosomes(df):
    """Filter out all variants that do not reside on the primary chromosomes
    (i.e. 1-22, X, Y, MT)."""

    # Print excluded variants for debugging/logging purpose.
    exclude_indices = df.index[~df['Chromosome'].isin(["%s" % chrom for chrom in range(1, 23)] + ['X', 'Y', 'MT'])]
    if len(exclude_indices) > 0:
        print "##"
        print "## WARNING: Filtering out all variants that do not reside on the primary chromosomes:"
        print "##"
        print df[df.index.isin(exclude_indices)]

    # Drop the variants
    df.drop(exclude_indices, inplace=True)

    return df

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_genomics", action="store", required=True,
                        metavar='FILE',
                        help="Path to input TCGA Genomics Data.")

    parser.add_argument("--out_maf", action="store", required=True,
                        metavar='FILE',
                         help="Path to the output MAF file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_genomics: ", repr(options.in_genomics)
    print "##   out_maf:", repr(options.out_maf)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
