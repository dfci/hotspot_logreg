#!/usr/bin/env python2

"""
SYNOPSIS
    Script to validate the Individual Cancer Centers' Target BED files.

    Use the BED specification defined at 'https://genome.ucsc.edu/FAQ/FAQformat#format1'
    as starting point. However, also add additional (more strict) validation
    conditions.

NOTES

    Individual BED file should be validated against this script before being
    uploaded/registered on IntelCCC's cluster.

    If Validation SUCCEEDED, script will exit with return code = 0.

    If Validation FAILED, script will exit with return code = 99.

    If an unexpected ERROR occur, the script should exit with other non-zero
    return code.

EXAMPLES

    # Validate the 'METABRIC' BED file.
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/public/metabric.targetedIntervals.r1.bed

    # Validate the 'DFCI-ONCOPANEL-1' BED file..
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/private/dfci.PROFILE_POPv1.filtered.r1.bed

    # Validate the 'DFCI-ONCOPANEL-2' BED file.
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/private/dfci.PROFILE_POPv2.filtered.r1.bed

    # Validate the 'DFCI-ONCOPANEL-3' BED file.
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/private/dfci.PROFILE_POPv3.r1.bed

    # Validate the 'GENIE-MSK-IMPACT341' BED file..
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/public/genie/genie_msk.IMPACT341.r1.bed

    # Validate the 'GENIE-MSK-IMPACT410' BED file.
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/public/genie/genie_msk.IMPACT410.r1.bed

    # Validate the 'WGS' BED file.
    #     Expected Exit Code: 0
    ./validate_bed.py --in_bed ../examples/inputs/public/wgs.target.r1.bed

EDGE_CASES

    # Edge Case #1: Run validation on Broken BED file that is purposely designed
    #               to exercise the code for syntax errors and ability to detect
    #               errors/warnings.
    #     Gracefully Exit with Return Code: 99
    ./validate_bed.py --in_bed ../examples/inputs/public/broken.bed; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #2: Input MAF file doesn't exist.
    #     Gracefully Exit with Return Code: 99
    ./validate_bed.py --in_bed ../examples/inputs/public/missing.bed; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #3: Pass in a Raw (Bed-like) Interval file. Test if code can gracefully handle 
    #               weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_bed.py --in_bed ../examples/inputs/public/raw/metabric.targetedIntervals.txt; \
    echo -e "##\n## RETURN_CODE: $?"

    # Edge Case #4: Pass in a MAF file. Test if code can gracefully handle weird input formats.
    #     Gracefully Exit with Return Code: 99
    ./validate_bed.py --in_bed ../examples/inputs/public/metabric.SNV.r4.maf; \
    echo -e "##\n## RETURN_CODE: $?"

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (Dec-2016)
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

VALIDATION_FAIL_RETURN_CODE = 99

def main(options):

    if not os.path.isfile(options.in_bed):
        sys.stderr.write("##\n")
        sys.stderr.write("## ERROR: BED file '%s' is not a valid existing file.\n" % options.in_bed)
        sys.exit(VALIDATION_FAIL_RETURN_CODE)

    err_toofewcols_list = list()
    err_chromname_list = list()
    err_nullvalue_list = list()
    err_nonint_chromstart_list = list()
    err_nonint_chromend_list = list()
    err_outofrange_chromstart_list = list()
    err_outofrange_chromend_list = list()
    err_diffpos_list = list()

    line_num = 0

    for line in open(options.in_bed, 'rU'):

        line_num += 1

        # Ignore all comment lines (i.e. line that starts with '#' character).
        if line.startswith("#"):
            continue

        # Only '\n' since opened file with universal newline support.
        line = line.rstrip('\n')

        # Assume the data is tab-delimited.
        cols = line.split('\t')

        if len(cols) < 3:
            err_toofewcols_list.append("Line #%d" % line_num)
            continue

        has_null_value = False
        for col_index in range(3):
            if cols[col_index].strip() == "":
                has_null_value = True

        if has_null_value:
            err_nullvalue_list.append("Line #%d" % line_num)
            continue

        # Assume that there are no column header line. Will just reference the
        # columns by index.
        chrom_name = cols[0]
        if not check_chrom_name(chrom_name):
            err_chromname_list.append("Line #%d: chromName->%s" % (line_num, repr(cols[0])))
            continue

        try:
            chrom_start = int(cols[1])
        except:
            err_nonint_chromstart_list.append("Line #%d: chromStart->%s" % (line_num, repr(cols[1])))
            continue

        try:
            chrom_end = int(cols[2])
        except:
            err_nonint_chromend_list.append("Line #%d: chromEnd->%s" % (line_num, repr(cols[2])))
            continue

        if not check_genomics_pos(chrom_name, chrom_start, offset=1):
            err_outofrange_chromstart_list.append("Line #%d: chromName->%s; chromStart->%s" % (line_num, repr(chrom_name), repr(chrom_start)))

        if not check_genomics_pos(chrom_name, chrom_end, offset=0):
            err_outofrange_chromend_list.append("Line #%d: chromName->%s; chromEnd>%s" % (line_num, repr(chrom_name), repr(chrom_end)))

        # Because chrom_start is 0-based and chrom_end is 1-based, use '>=' rather than '>'.
        if chrom_start >= chrom_end:
            err_diffpos_list.append("Line #%d: chromStart->%s; chromEnd->%s" % (line_num, repr(chrom_start), repr(chrom_end)))

    # Add error messages to ERRORS
    max_show = 2

    add_truncated_list_err_msg("less than 3 columns",
                               err_toofewcols_list, max_show)

    add_truncated_list_err_msg("null/missing value(s) ",
                               err_nullvalue_list, max_show)

    add_truncated_list_err_msg("invalid 'ChromName (1st col)' value",
                               err_chromname_list, max_show)

    add_truncated_list_err_msg("invalid 'ChromStart (2nd col)' value",
                               err_nonint_chromstart_list, max_show)

    add_truncated_list_err_msg("invalid 'ChromEnd (3rd col)' value",
                               err_nonint_chromend_list, max_show)

    add_truncated_list_err_msg("out-of-range 'ChromStart (2nd col)' value",
                               err_outofrange_chromstart_list, max_show)

    add_truncated_list_err_msg("out-of-range 'ChromEnd (3rd col)' value",
                               err_outofrange_chromend_list, max_show)

    add_truncated_list_err_msg("ChromStart >= ChromEnd",
                               err_diffpos_list, max_show)

    summarize_and_exit()

def check_chrom_name(chrom_name):
    """Check validity of the input chrom_name."""

    expected_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
                       '12', '13', '14', '15', '16', '17', '18', '19', '20',
                       '21', '22', 'X', 'Y', 'MT']

    if chrom_name in expected_chroms:
        return True
    else:
        return False

def check_genomics_pos(chrom_name, position, offset):
    """Check genomic position (either chrom_start or chrom_end)

    Note
    ----
    Set offset to 1 for chrom_start since, chrom_start is 0-based.
    Set offset to 0 for chrom_end since, chrom_end is 1-based.
    """
    
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

    if chrom_name not in chromosome_length:
        # err_msg already produced by check_chrom_name
        return True

    if position + offset <= 0:
        return False

    if position + offset > chromosome_length[chrom_name]:
        return False

    return True

def add_truncated_list_err_msg(msg_body, msg_list, max_show):
    """Create error_message string showing only the first 'max_show' error cases and
    add it to ERRORS global variable"""

    if len(msg_list) > 0:
        msg = "Found %s lines with %s. E.g. %s" % (len(msg_list), msg_body, msg_list[:max_show])
        ERRORS.append(msg)

def summarize_and_exit():
    """Summarize the BED validation result and exit

    Notes
    -----
    (1) If encountered no error, exit with return code = 0.

    (2) If encountered error(s), exit with return code = 99.
    """

    if len(ERRORS) != 0:
        sys.stderr.write("##\n")
        sys.stderr.write("## BED Validation FAILED. Please address the following critical ERROR(s):\n")
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
        sys.stderr.write("## BED Validation SUCCEEDED with the following WARNING(s):\n")
        for index, msg in enumerate(WARNINGS, start=1):
            sys.stderr.write("##\n")
            sys.stderr.write("##    %d. %s\n" % (index, msg))

        sys.exit(0)

    else:
        sys.stderr.write("##\n")
        sys.stderr.write("## BED Validation SUCCEEDED with no WARNING.\n")
        sys.exit(0)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_bed", action="store", required=True,
                        metavar='FILE',
                        help="Path to the input BED file for validation.")

    options = parser.parse_args()

    print "## Program Version:", repr(VERSION)
    print "## Validating BED file %s." % repr(options.in_bed)

    main(options)

    sys.exit(0)
