#!/usr/bin/env python2

"""
SYNOPSIS
    Script to validate the Individual Cancer Center MAF files.

NOTES

    Individual MAF file should be validated against this script before being
    uploaded/registered on IntelCCC's cluster.

    If Validation SUCCEEDED, script will exit with return code = 0.

    If Validation FAILED, script will exit with return code = 99.

    If an unexpected ERROR occur, the script should exit with other non-zero
    return code.

EXAMPLES

    # Validate the 'METABRIC' Full MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/public/metabric/metabric.maf

    # Validate the 'METABRIC' (Subset_Genes) MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/public/metabric/metabric.subset_genes.maf

    # Validate the 'DFCI-ONCOPANEL-1' MAF file..
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.maf

    # Validate the 'DFCI-ONCOPANEL-2' MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/private/dfci.DFCI-ONCOPANEL-2.r4.maf

    # Validate the 'DFCI-ONCOPANEL-3' MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/private/dfci.DFCI-ONCOPANEL-3.r4.maf

    # Validate the 'MSK-IMPACT341' MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/public/genie/genie_msk.IMPACT341.r3.maf

    # Validate the 'MSK-IMPACT410' MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/public/genie/genie_msk.IMPACT341.r3.maf

    # Validate the 'TCGA-BRCA' MAF file.
    #     Expected Exit Code: 0
    ./validate_maf.py --in_maf ../examples/inputs/public/tcga_brca/tcga.cell_2015.r2.maf

EDGE_CASES

    # Edge Case #1: Exercise code for syntax errors and ability to detect errors/warnings.
    #     Gracefully Exit with Return Code: 99
    ./validate_maf.py --in_maf ../examples/inputs/private/dfci.DFCI-ONCOPANEL-1.r4.maf --test_code; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #2: Input MAF file doesn't exist.
    #     Gracefully Exit with Return Code: 99
    ./validate_maf.py --in_maf ../examples/inputs/public/missing.maf; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #3: Pass in a Raw Variant Call file. Test if code can gracefully handle weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_maf.py --in_maf ../examples/inputs/public/raw/metabric.somaticMutations_incNC.txt; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #4: Pass in a BED file. Test if code can gracefully handle weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_maf.py --in_maf ../examples/inputs/public/metabric.targetedIntervals.r1.bed; \
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
import numpy as np
import pandas as pd

VERSION = 'mc_v0.9_dev'

ERRORS = list()
WARNINGS = list()

TEST_ROW = 0 # Use ONLY with --test_code option.
VALIDATION_FAIL_RETURN_CODE = 99

def main(options):

    maf_df = import_maf(options.in_maf)

    check_column_header(maf_df)

    check_missing_data(maf_df)

    check_chrom_values(maf_df)

    check_genomics_pos(maf_df)

    check_strand(maf_df)

    check_ncbi_build(maf_df)

    check_center(maf_df)

    check_tumor_seq_allele1(maf_df)

    check_alleles_and_positions(maf_df)

    summarize_and_exit()

def check_tumor_seq_allele1(maf_df):
    """Check that the 'Tumor_Seq_Allele1' value match 'Reference_Allele' for
    all variant calls."""

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], 'Reference_Allele'] = 'AAAACCCCGGGG'; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'Tumor_Seq_Allele1'] = 'AAAACCCCGGGG'; TEST_ROW+=1

    if 'Reference_Allele' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    if 'Tumor_Seq_Allele1' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    mismatch_count = len(maf_df[maf_df['Reference_Allele'] != maf_df['Tumor_Seq_Allele1']])

    if mismatch_count > 0:
        err_msg = "Mismatch between 'Reference_Allele' and 'Tumor_Seq_Allele1' in %d rows." % mismatch_count
        ERRORS.append(err_msg)

def check_chrom_values(maf_df):
    """Check that values in the 'Chromosome' column."""

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], 'Chromosome'] = 'chr1'; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'Chromosome'] = 'foochr'; TEST_ROW+=1

    expected_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
                       '12', '13', '14', '15', '16', '17', '18', '19', '20',
                       '21', '22', 'X', 'Y', 'MT']

    if 'Chromosome' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    chroms = list(maf_df['Chromosome'].unique())

    invalid_chroms = set(chroms) - set(expected_chroms)

    if len(invalid_chroms) > 0:
        err_msg = "Invalid value(s) in 'Chromosome' column: %s" % list(invalid_chroms)
        ERRORS.append(err_msg)

def check_genomics_pos(maf_df):
    """Check that 'Start_Position' and 'End_Position' values are in range for
    the given GRCh37's Chromosome."""

    if options.test_code:
        global TEST_ROW
        for pos_key in ['Start_Position', 'End_Position']:
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = "-100000"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = "1000000000"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = "-98765"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = "999999999"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = 'abc'; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], pos_key] = 'xyz'; TEST_ROW+=1

    if 'Chromosome' not in maf_df:
        # err_msg already produced by check_column_header().
        return

    # Derive chromosome_length from 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai'
    chromosome_length = dict()
    chromosome_length["1"] = 249250621
    chromosome_length["2"] = 243199373
    chromosome_length["3"] = 198022430
    chromosome_length["4"] = 191154276
    chromosome_length["5"] = 180915260
    chromosome_length["6"] = 171115067
    chromosome_length["7"] = 159138663
    chromosome_length["8"] = 146364022
    chromosome_length["9"] = 141213431
    chromosome_length["10"] = 135534747
    chromosome_length["11"] = 135006516
    chromosome_length["12"] = 133851895
    chromosome_length["13"] = 115169878
    chromosome_length["14"] = 107349540
    chromosome_length["15"] = 102531392
    chromosome_length["16"] = 90354753
    chromosome_length["17"] = 81195210
    chromosome_length["18"] = 78077248
    chromosome_length["19"] = 59128983
    chromosome_length["20"] = 63025520
    chromosome_length["21"] = 48129895
    chromosome_length["22"] = 51304566
    chromosome_length["X"] = 155270560
    chromosome_length["Y"] = 59373566
    chromosome_length["MT"] = 16569

    for pos_key in ['Start_Position', 'End_Position']:

        if pos_key not in list(maf_df):
            # err_msg already produced by check_column_header().
            continue

        nonint_pos_list = list()
        outofrange_pos_list = list()

        for index, row in maf_df.iterrows():

            if str(row[pos_key]) == 'nan':
                # err_msg already produced by check_missing_data().
                continue

            if str(row['Chromosome']) == 'nan':
                # err_msg already produced by check_missing_data().
                continue

            chrom = row['Chromosome']
            position = row[pos_key]

            try:
                position = int(position)
            except:
                nonint_pos_list.append("%s of Chrom %s" % (position, chrom))
                continue

            if position <= 0:
                outofrange_pos_list.append("%s of Chrom %s" % (position, chrom))

            if chrom in chromosome_length:
                if position > chromosome_length[chrom]:
                    outofrange_pos_list.append("%s of Chrom %s" % (position, chrom))

        if len(nonint_pos_list) > 0:
            err_msg = "Non-integer '%s' values: %s" % (pos_key, nonint_pos_list[:3])
            if len(nonint_pos_list) > 3:
                err_msg.rstrip("]")
                err_msg += "(MORE...)]"
            ERRORS.append(err_msg)

        if len(outofrange_pos_list) > 0:
            err_msg = "Out-of-range '%s' values: %s" % (pos_key, outofrange_pos_list[:3])
            if len(outofrange_pos_list) > 3:
                err_msg.rstrip("]")
                err_msg += "(MORE...)]"
            ERRORS.append(err_msg)

def check_alleles_and_positions(maf_df):
    """Check for consistencies between 'Ref_Allele', 'Tumor_Seq_Allele2',
    'Start_Position', and 'End_Position'."""

    if options.test_code:
        global TEST_ROW

        maf_df.loc[maf_df.index[TEST_ROW], 'Reference_Allele'] = "a";
        maf_df.loc[maf_df.index[TEST_ROW], 'Tumor_Seq_Allele2'] = "g"; TEST_ROW+=1

        maf_df.loc[maf_df.index[TEST_ROW], 'Start_Position'] = "100"
        maf_df.loc[maf_df.index[TEST_ROW], 'End_Position'] = "99"; TEST_ROW+=1

        for allele_key in ['Reference_Allele', 'Tumor_Seq_Allele2']:
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = ""; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = "g"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = "Cc"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = "hg19"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = "foo"; TEST_ROW+=1
            maf_df.loc[maf_df.index[TEST_ROW], allele_key] = "-"; TEST_ROW+=1

        maf_df.loc[maf_df.index[TEST_ROW], 'Reference_Allele'] = "-";
        maf_df.loc[maf_df.index[TEST_ROW], 'Tumor_Seq_Allele2'] = "cgt"; TEST_ROW+=1

        maf_df.loc[maf_df.index[TEST_ROW], 'Reference_Allele'] = "A";
        maf_df.loc[maf_df.index[TEST_ROW], 'Tumor_Seq_Allele2'] = "GG"; TEST_ROW+=1

    for colname in ['Reference_Allele', 'Tumor_Seq_Allele2', 'Start_Position', 'End_Position']:
        if colname not in maf_df:
            # err_msg already produced by check_column_header().
            return

    err_zerolenrefallele_list = list()
    err_zerolenaltallele_list = list()
    err_refallele_list = list()
    err_altallele_list = list()
    err_refaltalleles_list = list()

    err_negdiffpos_list = list()
    err_snpdiffpos_list = list()
    err_insdiffpos_list = list()
    err_deldiffpos_list = list()

    warn_complexvariant_list = list()

    for index, row in maf_df.iterrows():

        # Check Allele Consistency
        ref_allele = row['Reference_Allele']
        alt_allele = row['Tumor_Seq_Allele2']

        if str(ref_allele) == 'nan' or str(alt_allele) == 'nan':
            # err_msg already produced by check_missing_data().
            continue

        if len(ref_allele) == 0:
            err_zerolenrefallele_list.append("Row #%d: Ref_Alelle->%s" % (index, repr(ref_allele)))
            continue

        if len(alt_allele) == 0:
            err_zerolenaltallele_list.append("Row #%d: Tumor_Seq_Alelle2->%s" % (index, repr(alt_allele)))
            continue

        ref_allele_nonACGTchars = set(ref_allele) - set(["A", "C", "G", "T"])
        alt_allele_nonACGTchars = set(alt_allele) - set(["A", "C", "G", "T"])

        if ref_allele != "-" and len(ref_allele_nonACGTchars) > 0:
            err_refallele_list.append("Row #%d: Ref_Allele->%s" % (index, repr(ref_allele)))
            continue

        if alt_allele != "-" and len(alt_allele_nonACGTchars) > 0:
            err_altallele_list.append("Row #%d: Tumor_Seq_Allele2->%s" % (index, repr(alt_allele)))
            continue

        # Check Position + Allele Consistency
        try:
            start_pos = int(row['Start_Position'])
            end_pos = int(row['End_Position'])
        except:
            # err_msg already produced by check_genomics_pos()
            continue

        diff_pos = end_pos - start_pos
        if diff_pos < 0:
            err_negdiffpos_list.append("Row #%d: Start->%d; End->%d" % (index, start_pos, end_pos))
            continue

        if ref_allele == "-" and alt_allele == "-":
            err_refaltalleles_list.append("Row #%d: Ref_Allele:%s; Tumor_Seq_Allele2->%s" % (index, repr(ref_allele), repr(alt_allele)))

        elif alt_allele == "-":
            # DELETION
            # e.g. ref_allele 'G' and alt_allele = '-'
            # e.g. ref_allele 'GCTT' and alt_allele = '-'
            if (diff_pos + 1) != len(ref_allele):
                err_deldiffpos_list.append("Row #%d: Start->%d; End->%d: Ref_Allele->%s" % (index, start_pos, end_pos, repr(ref_allele)))

        elif ref_allele == "-":
            # INSERTION
            # e.g. ref_allele '-' and alt_allele = 'G'
            # e.g. ref_allele '-' and alt_allele = 'GCTT'
            if diff_pos != 1:
                err_insdiffpos_list.append("Row #%d: Start->%d; End->%d" % (index, start_pos, end_pos))

        elif len(ref_allele) == len(alt_allele):
            # SNP, DNP, TNP, or ONP
            # e.g. ref_allele = 'A' and alt_allele = 'C'
            # e.g. ref_allele = 'AAAA' and  alt_allele = 'CCCC'
            if (diff_pos + 1) != len(ref_allele):
                err_snpdiffpos_list.append("Row #%d: Start->%d; End->%d: Ref_Allele->%s" % (index, start_pos, end_pos, repr(ref_allele)))

        else: # len(ref_allele) != len(alt_allele)
            # Potential a DELINs (combination of Deletion + Insertions)
            # e.g. ref_allele = 'T' and alt_allele = 'GCAAT'
            # e.g. ref_allele = 'GCAAT' and alt_allele = 'T'
            # Report them as WARNINGS and don't perform any check at the moment.
            warn_complexvariant_list.append("Row #%d: Start->%d; End->%d: Ref_Allele->%s Tumor_Seq_Allele2->%s" % (index, start_pos, end_pos,
                                                                                                                   repr(ref_allele),
                                                                                                                   repr(alt_allele)))

    # Add error messages to ERRORS
    add_truncated_list_err_msg("'len(Ref_Allele) == 0'",
                               err_zerolenrefallele_list, max_show=2)

    add_truncated_list_err_msg("'len(Tumor_Seq_Allele2) == 0'",
                               err_zerolenaltallele_list, max_show=2)

    add_truncated_list_err_msg("Invalid 'Ref_Allele' values",
                               err_refallele_list, max_show=2)

    add_truncated_list_err_msg("Invalid 'Tumor_Seq_Allele2' values",
                               err_altallele_list, max_show=2)

    add_truncated_list_err_msg("Invalid 'Ref_Allele' and 'Tumor_Seq_Allele2' combination",
                               err_refaltalleles_list, max_show=2)

    add_truncated_list_err_msg("'Start_Position > End_position'",
                               err_negdiffpos_list, max_show=2)

    add_truncated_list_err_msg("(SNP/MNP) Discrepancies between Start_Pos, End_Pos, and len(Ref_Allele)",
                               err_snpdiffpos_list, max_show=2)

    add_truncated_list_err_msg("(INSERTION) 'Start_Pos != End_Pos + 1'",
                               err_insdiffpos_list, max_show=2)

    add_truncated_list_err_msg("(DELETION) Discrepancies between Start_Pos, End_Pos, and len(Ref_Allele)",
                               err_deldiffpos_list, max_show=2)

    add_truncated_list_warn_msg("Found Complex Variants (e.g. Deletion+Insertion)",
                               warn_complexvariant_list, max_show=2)

def add_truncated_list_err_msg(msg_prefix, msg_list, max_show):
    """Create error_message string showing only the first 'max_show' error cases and
    add it to ERRORS global variable"""

    if len(msg_list) > 0:
        msg = msg_prefix
        msg += " in %s Rows. E.g. %s" % (len(msg_list), msg_list[:max_show])
        ERRORS.append(msg)

def add_truncated_list_warn_msg(msg_prefix, msg_list, max_show):
    """Create warning_message string showing only the first 'max_show' warning cases and
    add it to WARNINGS global variable"""

    if len(msg_list) > 0:
        msg = msg_prefix
        msg += " in %s Rows. E.g. %s" % (len(msg_list), msg_list[:max_show])
        WARNINGS.append(msg)

def check_strand(maf_df):
    """Check that the 'Strand' value of every row is '+'."""

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], 'Strand'] = "-"; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'Strand'] = "foo"; TEST_ROW+=1

    if 'Strand' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    valid_values = ["+"]
    values = list(maf_df['Strand'].unique())
    invalid_values = set(values) - set(valid_values)

    if len(invalid_values) > 0:
        err_msg = "Variant should always be reported on the '+' strand. Invalid values in 'Strand' column: %s" % list(invalid_values)
        ERRORS.append(err_msg)

def check_ncbi_build(maf_df):
    """Check that the 'NCBI_Build' value of each row is either 'hg19' or 'GRCh37'."""

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], 'NCBI_Build'] = "hg18"; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'NCBI_Build'] = "36.1"; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'NCBI_Build'] = "hg19"; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'NCBI_Build'] = "GRCh37"; TEST_ROW+=1

    if 'NCBI_Build' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    valid_values = ['hg19', 'GRCh37']
    values = list(maf_df['NCBI_Build'].unique())
    invalid_values = set(values) - set(valid_values)
    intersect_values = set(values) & set(valid_values)

    if len(invalid_values) > 0:
        err_msg = "Reference Assembly should be either 'hg19' or 'GRCh37'. Invalid values in 'NBCI_Build' column: %s" % list(invalid_values)
        ERRORS.append(err_msg)

    if len(intersect_values) > 1:
        warn_msg = "Found both 'hg19' and 'GRCh37' Reference Assembly ('NBCI_Build') values."
        WARNINGS.append(warn_msg)

def check_center(maf_df):
    """Check for consistency in the 'Center' value in the rows."""

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], 'Center'] = 'Foo'; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'Center'] = 'Bar'; TEST_ROW+=1
        maf_df.loc[maf_df.index[TEST_ROW], 'Center'] = ''; TEST_ROW+=1

    if 'Center' not in list(maf_df):
        # err_msg already produced by check_column_header().
        return

    values = list(maf_df['Center'].unique())

    if len(values) > 1:
        err_msg = "Expect single value in 'Center' column across all rows. Found multiple values in 'Center' column: %s" % list(values)
        ERRORS.append(err_msg)

def check_column_header(maf_df):
    """Check for missing or extra MAF columns in header."""

    if options.test_code:
        maf_df.drop('Score', axis=1, inplace=True)
        maf_df.drop('Matched_Norm_Sample_UUID', axis=1, inplace=True)

    expected_columns = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
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

    columns = list(maf_df)

    extra_columns = list(set(columns) - set(expected_columns))
    missing_columns = list(set(expected_columns) - set(columns))

    if len(extra_columns) > 0:
        err_msg = "Extra column(s) in MAF header: %s" % extra_columns
        ERRORS.append(err_msg)

    if len(missing_columns) > 0:
        err_msg = "Missing expected column(s) in MAF header: %s" % missing_columns
        ERRORS.append(err_msg)

def check_missing_data(maf_df):
    """Check for missing data in agreed 'non_null' columns."""

    non_null_columns = ["Chromosome", "Start_Position", "End_Position",
                        "Strand", "Reference_Allele", "Tumor_Seq_Allele1",
                        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Center",
                        "NCBI_Build"]

    if options.test_code:
        global TEST_ROW
        maf_df.loc[maf_df.index[TEST_ROW], "Start_Position"] = np.nan; TEST_ROW+=1

        for column in non_null_columns:
            maf_df.loc[maf_df.index[TEST_ROW], column] = np.nan
        TEST_ROW+=1

        for column in non_null_columns:
            maf_df.loc[maf_df.index[TEST_ROW], column] = np.nan; TEST_ROW+=1

        for column in list(maf_df):
            maf_df.loc[maf_df.index[TEST_ROW], column] = np.nan
        TEST_ROW+=1

    err_msg_list = list()

    for column in non_null_columns:

        if column not in list(maf_df):
            # err_msg already produced by check_column_header().
            continue

        null_counts =  maf_df[column].isnull().values.sum()
        if null_counts != 0:
            err_msg_list.append([column, null_counts])

    if len(err_msg_list) != 0:
        err_msg = "Missing data in column(s): ["
        err_msg += ", ".join(["%s(rows=%d)" % (repr(x[0]), x[1]) for x in err_msg_list])
        err_msg += "]"
        ERRORS.append(err_msg)

def import_maf(infile):
    """Import the MAF file as a pandas data frame

    Notes
    -----
    (1) This infile is a tab-delimited MAF file. For full-spec, see:

            https://github.com/mskcc/vcf2maf/blob/master/docs/vep_maf_readme.txt
    """

    try:
        maf_df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header = 0)
    except Exception as E:
        sys.stderr.write("##\n")
        sys.stderr.write("## ERROR: Fail to import MAF file: %s\n" % repr(E))
        sys.exit(VALIDATION_FAIL_RETURN_CODE)

    return maf_df

def summarize_and_exit():
    """Summarize the MAF validation result and exit (with status 1 if found
    errors)."""

    if len(ERRORS) != 0:
        sys.stderr.write("##\n")
        sys.stderr.write("## MAF Validation FAILED. Please address the following critical ERROR(s):\n")
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
        sys.stderr.write("## MAF Validation SUCCEEDED with the following WARNING(s):\n")
        for index, msg in enumerate(WARNINGS, start=1):
            sys.stderr.write("##\n")
            sys.stderr.write("##    %d. %s\n" % (index, msg))

        sys.exit(0)

    else:
        sys.stderr.write("##\n")
        sys.stderr.write("## MAF Validation SUCCEEDED with no WARNING.\n")
        sys.exit(0)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input MAF file for validation.")

    parser.add_argument("--test_code", action='store_true',
                        help="Exercise code for syntax errors and ability to detect errors/warnings.")

    options = parser.parse_args()

    print "## Program Version:", repr(VERSION)
    print "## Validating MAF file %s." % repr(options.in_maf)

    main(options)

    sys.exit(0)
