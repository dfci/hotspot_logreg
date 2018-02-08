#!/usr/bin/env python2

"""
SYNOPSIS
    Preprocess the METABRIC Breast Cancer data into a MAF format that is suitable
    for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Assume that "Reference_Allele"/"Tumor_Seq_Allele1"/"Tumor_Seq_Allele2"
        values in input ('--in_genomics') already satisfy the MAF-specification.

    (2) Import the 'details' preprocessed clinical data file. This clinical data
        contains:
            (A) A 'Exclude_Sample' column which indicate whether the sample
                would be excluded (filtered-out) from the outfile based on
                clinical filtering citeria.

                Note
                ----
                See 'prep_metabric_clinical.py' script for more information regarding
                the clinical filtering citeria to exclude samples.


EXAMPLES

    # (1) Default
    ./prep_metabric_maf.py \
        --in_genomics ../public/metabric/raw/brca_metabric/data_mutations_extended.txt \
        --in_clinical ../public/clinical/metabric.clinical_data.details.tsv \
        --out_maf ../public/metabric/metabric.maf

    # (2) Filter and output only variants on PIK3CA and AKT1 genes.
    ./prep_metabric_maf.py \
        --in_genomics ../public/metabric/raw/brca_metabric/data_mutations_extended.txt \
        --in_clinical ../public/clinical/metabric.clinical_data.details.tsv \
        --select_genes PIK3CA AKT1 \
        --out_maf ../public/metabric/metabric.subset_genes.maf

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

    maf_sequenced_samples_line = open(options.in_genomics, 'rU').readlines()[0]
    maf_sample_barcodes = maf_sequenced_samples_line.lstrip('#Sequenced_Samples:').rstrip('\n').split()

    genomics_df = import_metabric_genomics(options.in_genomics, options.select_genes)
    genomics_df = filter_non_primary_chromosomes(genomics_df)

    clinical_df = import_prepped_clinical(options.in_clinical)

    # Check that the maf_sample_barcodes exactly match samples in the clinical_df
    assert set(maf_sample_barcodes) == set(clinical_df['Tumor_Sample_Barcode'].unique())

    # Check that the samples in genomics_df is subset of maf_sample_barcodes
    assert set(genomics_df['Tumor_Sample_Barcode'].unique()) <= set(maf_sample_barcodes)

    # Filter-out 'Exclude_Sample' (see 'prep_metabric_clinical.py' for details).
    filtered_clinical_df = clinical_df[clinical_df['Exclude_Sample'] == 'False']

    # Merge the filtered clinical and genomics df.
    out_df = pd.merge(genomics_df,
                      filtered_clinical_df,
                      how="inner",  # INNER JOIN
                      left_on="Tumor_Sample_Barcode",
                      right_on="Tumor_Sample_Barcode",
                      sort=False,
                      indicator='indicator_column')


    # Check for duplicate rows
    any_duplicates = out_df.duplicated(subset=["Tumor_Sample_Barcode", "Start_Position",
                                               "End_Position", "Reference_Allele",
                                               "Tumor_Seq_Allele2"]).any()

    assert any_duplicates == False

    # Check that the values in these two columns are the same for every row.
    assert (out_df['Tumor_Seq_Allele1'] == out_df['Reference_Allele']).all()

    # Set Center to 'METABRIC'
    out_df['Center'] = 'METABRIC'

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
          (1)  Tumor_Sample_Barcode (e.g. 'MB-7115')
          (3)  Exclude_Sample ('True', 'False')
    """

    df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header=0)

    # Check that all the require columns exist
    required_columns = ['Tumor_Sample_Barcode', 'Exclude_Sample']

    assert set(required_columns) <= set(df)
    df = df[required_columns]

    # Ensure that there is no missing data in any of the columns.
    assert not df.isnull().values.any()

    # Ensure that there are no duplicated 'Tumor_Sample_Barcode' values
    assert not df.duplicated(subset=['Tumor_Sample_Barcode']).any()

    # Ensure that there is no unexpected 'Exclude_Sample' values.
    assert set(df['Exclude_Sample'].unique()) <= set(['True', 'False'])

    return df

def import_metabric_genomics(infile, select_genes):
    """Import the METABRIC genomics file (which actually is itself a MAF file
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

    if select_genes != None:
        # Keep only the rows that match specific genes
        df = df[df['Hugo_Symbol'].isin(select_genes)]

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

    # Ensure that there is no unexpected value of NCBI_Build column.
    assert df["NCBI_Build"].unique() == ["GRCh37"]

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
                        help="Path to input METABRIC Genomics Data.")

    parser.add_argument("--in_clinical", action="store", required=True,
                        metavar='FILE',
                        help="Path to input (Detailed Proprocessed) METABRIC Clinical Data.")

    parser.add_argument("--select_genes", action="store", default=None,
                        nargs='+',
                        help="Filter and output data for only specific genes")

    parser.add_argument("--out_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output MAF file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_genomics:", repr(options.in_genomics)
    print "##   in_clinical: ", repr(options.in_clinical)
    print "##   select_genes:", repr(options.select_genes)
    print "##   out_maf:", repr(options.out_maf)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
