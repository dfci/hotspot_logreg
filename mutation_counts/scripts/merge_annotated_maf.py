#!/usr/bin/env python2

"""
SYNOPSIS
    Merge data from multiple input annotated MAF files into a single output
    annotated MAF file.

NOTES

    (a) Assume that every input MAF file is annotated with maf2maf/vcf2maf v1.6.9
        (VEP 85) and has the following 113 columns:

            https://github.com/mskcc/vcf2maf/blob/v1.6.9/docs/vep_maf_readme.txt

        To simplify the merging logic, be a little strict here and raise an error
        if the input MAF is missing any of these 113 columns. The strict check
        also indirectly help ensure that the input MAF was not annotated with a
        different version of maf2maf/VEP.

        Each input MAF file can have additional columns (to these 113), but the
        additional columns will be discarded (i.e. not propagated to the output MAF).

    (b) Add a non-standard column ('Source_Maf') to the output MAF file to
        retain information about the source of each sample/variant_call:

        - Use the basename of the source input_maf file as the 'SEQ_ASSAY_ID'
          value of each row.

        - Also check that there is no collision between basenames of the
          different input_maf files.

        - If the input MAF file already have a 'Source_Maf' column, then the
          existing data in that column will be discarded/overridden.

    (c) Append 'Center' name to 'Tumor_Sample_Barcode' as a way to help prevent
        samples from different input MAF files from colliding. Motivation is to
        ensure that 'Tumor_Sample_Barcode' can be used on a unique identifier of
        samples in downstream analysis steps.

        - After modifing the 'Tumor_Sample_Barcode' values, then check that
          there is no collision between the 'Tumor_Sample_Barcode' values of
          the different input MAF files.

EXAMPLES

    # Merge 4 input MAF files (three from DFCI and one from METABRIC).
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                  ../examples/inputs/private/dfci.DFCI-ONCOPANEL-2.r4.vep.maf \
                  ../examples/inputs/private/dfci.DFCI-ONCOPANEL-3.r4.vep.maf \
                  ../examples/inputs/public/metabric/metabric.vep.maf \
        --out_maf private/dfci_plus_metabric.r6.vep.maf

    # Merge DFCI, MSK-IMPACT, METABRIC, TCGA, and SangerWGS datasets.
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                  ../examples/inputs/private/dfci.DFCI-ONCOPANEL-2.r4.vep.maf \
                  ../examples/inputs/private/dfci.DFCI-ONCOPANEL-3.r4.vep.maf \
                  ../examples/inputs/public/genie/genie_msk.IMPACT341.r3.vep.maf \
                  ../examples/inputs/public/genie/genie_msk.IMPACT410.r3.vep.maf \
                  ../examples/inputs/public/metabric/metabric.vep.maf \
                  ../examples/inputs/public/from_intelccc/sanger_somatic_coding_annotated.maf \
                  ../examples/inputs/public/tcga_brca/tcga.cell_2015.r2.vep.maf \
        --out_maf private/merged.dfci_plus_public.r6.vep.maf

    # Merge data from only Public sources (MSK-IMPACT, METABRIC, TCGA, and SangerWGS).
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/public/genie/genie_msk.IMPACT341.r3.vep.maf \
                  ../examples/inputs/public/genie/genie_msk.IMPACT410.r3.vep.maf \
                  ../examples/inputs/public/metabric/metabric.vep.maf \
                  ../examples/inputs/public/from_intelccc/sanger_somatic_coding_annotated.maf \
                  ../examples/inputs/public/tcga_brca/tcga.cell_2015.r2.vep.maf \
        --out_maf private/merged.only_public.r6.vep.maf

EDGE_CASES

    Edge Case #1: Single input MAF file.
        Expected Behavior: Gracefully handle and generate correct output.
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
        --out_maf private/tmp.merged.vep.maf

    Edge Case #1: Some of the input MAF files does not exist.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                  ../examples/inputs/public/missing.vep.maf \
        --out_maf private/tmp.merged.vep.maf

    Edge Case #2: Collision between basenames of the different input_maf files.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                  ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
        --out_maf private/tmp.merged.vep.maf

    Edge Case #3: Collision between the 'Tumor_Sample_Barcode' values of the different input MAF files.
        Expected Behavior: Raise Appropriate Exception.
    ./merge_annotated_maf.py \
        --in_mafs ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                  ../examples/inputs/private/dupl.dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
        --out_maf private/tmp.merged.vep.maf

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

# Expected columns of a MAF file annotated with maf2maf/vcf2maf v1.6.9 (VEP 85).
# See: https://github.com/mskcc/vcf2maf/blob/v1.6.9/docs/vep_maf_readme.txt
VEP_MAF_columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position',
                   'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele',
                   'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
                   'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
                   'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
                   'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status',
                   'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer',
                   'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID',
                   'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count',
                   'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position',
                   'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM',
                   'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS',
                   'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON',
                   'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF',
                   'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS',
                   'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO',
                   'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN',
                   'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'flanking_bps',
                   'variant_id', 'variant_qual', 'ExAC_AF_Adj']

def main(options):

    in_maf_dfs = list()

    in_maf_basenames = list()  # For collision check
    in_maf2samples = dict()  # For collision check

    for in_maf_file in options.in_mafs:

        in_maf_basename = os.path.basename(in_maf_file)
        in_maf_df = import_maf(in_maf_file)

        check_column_header(in_maf_df, in_maf_file)

        # Ensure no collision between basenames of the different input MAF files.
        if in_maf_basename in in_maf_basenames:
            raise Exception("Two input_maf files have the same basename: %s" % repr(in_maf_basename))
        in_maf_basenames.append(in_maf_basename)

        # Add a non-standard column ('Source_Maf') to the output MAF file to
        # retain information about the source of each sample/variant_call.
        in_maf_df['Source_Maf'] = in_maf_basename

        # Append 'Center' name to 'Tumor_Sample_Barcode' as a way to help prevent
        # sample_names from different centers from colliding.
        assert not in_maf_df["Center"].isnull().values.any()
        assert not in_maf_df["Tumor_Sample_Barcode"].isnull().values.any()
        in_maf_df["Tumor_Sample_Barcode"] = in_maf_df["Center"] + "_" + in_maf_df["Tumor_Sample_Barcode"]

        check_tumor_sample_barcode_collisions(in_maf_df, in_maf_file, in_maf2samples)

        # Reorder and retain only the VEP_MAF_columns + 'Source_Maf'.
        in_maf_df = in_maf_df[VEP_MAF_columns + ['Source_Maf']]

        in_maf_dfs.append(in_maf_df)

    # Concatenate and output the MAF.
    out_df = pd.concat(in_maf_dfs, axis=0)
    fout = open(options.out_maf, 'w')
    fout.write("#version 2.4\n")
    out_df.to_csv(fout, sep="\t", na_rep='', index=False)
    fout.close()

    # Print Merged MAF Summary for debugging/logging purposes.
    print "##", "-" * 50
    print "## Merged MAF Summary:"
    print "##   Total # Samples:", format(len(out_df['Tumor_Sample_Barcode'].unique()), ',d')
    print "##   Total # Variant Calls:", format(len(out_df), ',d')
    print "##", "-" * 50

def check_tumor_sample_barcode_collisions(curr_maf_df, curr_maf_file, prev_maf2samples):
    """Ensure no collision between 'Tumor_Sample_Barcode' values from 'curr_maf_df'
    and the previously processed in_maf_files ('prev_maf2samples').

    Notes
    -----
    (1) If there is collision, this function will raise an Exception.

    (2) If there is no collision, this function will add the samples from
        'curr_maf_df' to 'prev_maf2samples' dict and then return.
    """

    samples = set(curr_maf_df["Tumor_Sample_Barcode"].unique())
    assert curr_maf_file not in prev_maf2samples

    for prev_maf_file, prev_samples in prev_maf2samples.iteritems():
        duplicates = list(prev_samples.intersection(samples))

        if len(duplicates) != 0:
            max_show = max(len(duplicates), 5)
            err_msg = "Found %s duplicated 'Tumor_Sample_Barcode' values between" % len(duplicates)
            err_msg += " maf_1 %s and maf_2 %s." % (repr(curr_maf_file), repr(prev_maf_file))
            if len(duplicates) > 5:
                err_msg += " Here is a subset of the duplicated values: %s" % duplicates[:5]
            else:
                err_msg += " Here are the duplicated values: %s" % duplicates

            raise Exception(err_msg)

    prev_maf2samples[curr_maf_file] = samples

def check_column_header(maf_df, in_maf_file):
    """Check that in_maf_file have all 34 standard MAF columns."""

    columns = list(maf_df)

    missing_columns = list(set(VEP_MAF_columns) - set(columns))

    if len(missing_columns) > 0:
        raise Exception("Missing VEP/maf2maf MAF column(s) %s in input_maf %s." % (missing_columns,
                                                                                   repr(in_maf_file)))
def import_maf(infile):
    """Import the MAF file as a pandas data frame

    Notes
    -----
    (1) This infile is a tab-delimited MAF file. For full-spec, see:

            https://github.com/mskcc/vcf2maf/blob/master/docs/vep_maf_readme.txt
    """

    maf_df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header = 0)

    return maf_df

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_mafs", action="store", required=True, nargs='+',
                        metavar='FILE',
                        help="List of paths to input MAF files.")

    parser.add_argument("--out_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to the output (merged) MAF file.")

    options = parser.parse_args()

    for index in range(len(options.in_mafs)):
        options.in_mafs[index] = os.path.abspath(options.in_mafs[index])

    options.out_maf = os.path.abspath(options.out_maf)

    print "##", "-" * 50
    print "## Program Version:", repr(VERSION)
    print "## Specified Options:"
    print "##   in_mafs:"
    for index, in_maf in enumerate(options.in_mafs, start=1):
        print "##       (%d) %s" % (index, repr(in_maf))
    print "##   out_maf:", repr(options.out_maf)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
