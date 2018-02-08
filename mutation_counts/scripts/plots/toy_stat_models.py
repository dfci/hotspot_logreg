#!/usr/bin/env python2

"""
SYNOPSIS

    Run various toy (simplified) statistical tests and compare the resulting
    p_values generated over a range of PIK3CA_G118 mutation count values.

EXAMPLE:

    # (1) Run on the Toy Dataset
    ./toy_stat_models.py \
        --dataset_name TOY \
        --out_table p_values.TOY.PIK3CA_G118.txt \
        --out_plot plot.TOY.PIK3CA_G118.png

    # (2) Run on the slightly more slight realistics counts from METABRIC dataset
    ./toy_stat_models.py \
        --dataset_name METABRIC \
        --out_table p_values.METABRIC.PIK3CA_G118.txt \
        --out_plot plot.METABRIC.PIK3CA_G118.png

NOTES

    (1) The Toy Statistic Tests Performed are:

        (A) Pearson Chi-Square Test (with Correction)
        (B) Fisher's Exact Test
        (C) Univariate Logistic Regression (Likelihood Ratio Test) [Note: No Fixed Effects]
        (D) Univariate Logistic REgression (Wald's Test)           [Note: No Fixed Effects]

    (2) Toy Dataset:

        Consider a hypothetical situation where 5000 Bresat Cancer Samples were sequenced
        using a 10 genes panel. For simplicity, assume that each gene has a coding region
        of exactly 200 codon positions.

        Now suppose that PIK3CA is one of the gene in the panel. Further suppose that
        we observed 15 (out of 500 samples) with mutation call at codon G118 of PIK3CA gene.

        Question: Is PIK3CA_G118 a hotpost

        Question: Is PIK3CA_G118 mutated more frequently than the back rate.

        For simplicity, further assume that all other codons in the 10 gene panel
        (aside from PIK3CA_G118) mutates at the common background rate of 1 mutation per 100
        codons for all samples.

            - No variability in the background mutation rate across genes
            - No variability in the background mutation rate across codon types
            - No variability in the background mutation rate across samples.

        A = # of mutation calls observed at PIK3CA_G118 codon (across all samples)
        A = 15

        B = # of reference calls observed at PIK3CA_G118 codon (across all samples)
        B = 500 - 15
        B = 485

        C = # of mutation calls observed at all other codons in panel (across all samples)
        C = (# codon positiosn in panel - 1) * (background mutation rate) * (total # of samples)
        C = (10 genes * 200 codons/gene - 1) * (1 mutations per 100 codons) * 500
        C = 9,995

        D = # of reference calls observed at all other codons in panel (across all samples)
        D = (# codon positiosn in panel - 1) * (1.0 - background mutation rate) * (total # of samples)
        D = (10 genes * 200 codons/gene - 1) * (0.99) * 500
        D = 989,505

        Here is corresponding 2x2 contigency table.

                                                            Response Variable

                                                    Has Mutation    No Mutation
                                                  -------------------------------
                                                  |              |              |
                             Is At PIK3CA_G118    |    A = 15    |   B = 485    |
        Explanatory Variable                      |--------------|--------------|
                             Not At PIK3CA_G118   |              |              |
                                                  |   C = 9995   |  D = 989,505 |
                                                  |-----------------------------|


    (3) Slightly More Realistic counts observed in METABRIC dataset:

                                                            Response Variable

                                                    Has Mutation    No Mutation
                                                  ------------------------------------
                                                  |              |                   |
                             Is At PIK3CA_G118    |  A = 10      |  B = 6,846        |
        Explanatory Variable                      |--------------|-------------------|
                             Not At PIK3CA_G118   |              |                   |
                                                  | C = 10,096   | D = 1,734,285,348 |
                                                  |----------------------------------|

"""

import sys
import os
import time
import argparse
import pandas as pd
import numpy as np
import scipy as sp
import traceback
import random

import statsmodels.api as sm
import matplotlib.pyplot as plt

def main(options):

    if options.dataset_name == "TOY":
        A_plus_B = 500  # Num_Samples
        C = 9995
        num_genes = 10
        num_codons_per_gene = 200
        total_count = A_plus_B * num_genes * num_codons_per_gene

    elif options.dataset_name == "METABRIC":
        A_plus_B = 6846 + 10
        C = 10096
        total_count = 1734285348 + 10096 + 6846 + 10

    else:
        raise Exception("Invalid Dataset Name (%s)." % options.dataset_name)

    out_data_list = list()

    A_count_min = 0
    A_count_max = 20

    for A in range(A_count_min, A_count_max+1):

        print "##"
        print "## Set A value to %s." % A

        B = A_plus_B - A
        D = total_count - C - B - A

        if A == 0 or B == 0:
            continue

        attempt_num = 0
        while True:
            try:
                attempt_num += 1
                num_parts = attempt_num + 4
                p_value_pearson, p_value_fisher, p_value_lrt, p_value_wald = run_multiple_tests(A, B, C, D, num_parts)
            except:
                print "##"
                print "## WARNING: Logistic Regression failed to converge for attempt #%d. Here is the error message:" % attempt_num
                print "##"
                print traceback.print_exc()
                print "##"
                print "## WARNING: Retrying ..."
                print "##"
                if attempt_num == 10:
                    raise Exception("Logistic Regression failed to converge after %d attempts!" % attempt_num)
            else:
                if attempt_num > 1:
                    print "##"
                    print "## Successfully Ran Logistic Regression for after multiple attempts (#%d) (num_parts = %d)." % (attempt_num, num_parts)
                break

        out_data = dict()
        out_data["A_Count"] = A
        out_data["B_Count"] = B
        out_data["C_Count"] = C
        out_data["D_Count"] = D
        out_data["Pearson_Corrected"] = p_value_pearson
        out_data["Fisher"] = p_value_fisher
        out_data["Logistic_LRT"] = p_value_lrt
        out_data["Logistic_Wald"] = p_value_wald

        out_data_list.append(out_data)

    create_output_table(out_data_list, options.out_table)

    plot_results(options.out_table, options.dataset_name, options.out_plot)

def create_output_table(out_data_list, outfile):
    """Write the output to disk."""

    fout = open(options.out_table, 'w')

    header_cols = ["A_Count", "B_Count", "C_Count", "D_Count", "Pearson_Corrected", "Fisher", "Logistic_LRT", "Logistic_Wald"]

    fout.write("\t".join(header_cols) + "\n")

    for index in range(len(out_data_list)):

        out_data = out_data_list[index]

        cols = ["%d" % out_data["A_Count"],
                "%d" % out_data["B_Count"],
                "%d" % out_data["C_Count"],
                "%d" % out_data["D_Count"],
                "%E" % out_data["Pearson_Corrected"],
                "%E" % out_data["Fisher"],
                "%E" % out_data["Logistic_LRT"],
                "%E" % out_data["Logistic_Wald"]]

        line = "\t".join(cols)
        fout.write(line + '\n')

    fout.close()

def plot_results(in_table_path, dataset_name, outfile):
    """Import the data in the table and plot it."""

    my_df = pd.read_table(in_table_path, sep="\t", comment="#", header = 0)

    params = {'legend.fontsize': 10, 'axes.labelsize':13, 'axes.labelcolor':  "black",
              'axes.labelweight': "bold", "xtick.labelsize": 12, "ytick.labelsize": 12,
              'legend.markerscale': 1.0}
    plt.rcParams.update(params)

    if dataset_name == "TOY":
        x_min = 0
        x_max = 20.5
        y_min = 1E-11
        y_max = 1E+03
    elif dataset_name == "METABRIC":
        x_min = 0
        x_max = 20.5
        y_min = 1E-50
        y_max = 1E+02

    else:
        raise Exception("Invalid infile %s." % options.infile)

    point_size_fisher = 50
    print_size_LRT = 70
    point_size_other = 70

    plt.scatter(my_df['A_Count'], my_df['Fisher'],
                label="Fisher's Exact Test (2x2 Table)",
                s=point_size_fisher, color='green', zorder=10)

    plt.scatter(my_df['A_Count'], my_df['Pearson_Corrected'],
                label="Pearson's " + r'$\chi^2$' + " (2x2 Table)",
                s=point_size_other, color='blue', zorder=10)

    plt.scatter(my_df['A_Count'], my_df['Logistic_Wald'],
                label='Logistic Regression (two-sided Z-test on ' + r'$\beta_1$' +')',
                s=point_size_other, color='red', zorder=10)

    plt.scatter(my_df['A_Count'], my_df['Logistic_LRT'],
                label='Logistic Regression (likelihood-ratio test on ' + r'$\beta_1$' +')',
                s=print_size_LRT, color='magenta', zorder=5)

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.legend(loc="upper right")
    # plt.legend(loc="lower left")

    plt.title("Comparing PValues", fontsize = 20, fontweight='bold', y=1.02)

    plt.xlabel("Mutation Counts observed at PIK3CA_G118 (i.e. A Counts)")

    plt.ylabel("PValue", fontsize = 16)
    plt.yscale('log')
    plt.savefig(outfile)
    plt.clf()

def run_multiple_tests(A, B, C, D, num_parts):

    odds_1 = A / float(B)
    odds_2 = C / float(D)
    odds_ratio = odds_1 / odds_2

    print "## Counts:"
    print "##   A: %d" % A
    print "##   B: %d" % B
    print "##   C: %d" % C
    print "##   D: %d" % D

    print "##"
    print "## Odds:"
    print "##   odds_1: %f" % odds_1
    print "##   odds_2: %f" % odds_2
    print "##   odds_ratio: %f" % odds_ratio

    print "##"
    print "## Pearson 2x2 Chi-Square test:"
    # Pearson 2x2 Chi-square test
    obs = np.array([[A, B], [C, D]])
    chi2_pearson, _, dof, expected = sp.stats.chi2_contingency(obs, correction=True)

    assert dof == 1
    p_value_pearson = sp.stats.chi2.sf(chi2_pearson, dof)

    print "##   chi2 (Pearson):", repr(chi2_pearson)
    print "##   P-value (Pearson):", p_value_pearson
    print "##   dof:", repr(dof)
    print "##   expected:", repr(expected)

    print "##"
    print "## Fisher 2x2 Exact Test:"
    _,  p_value_fisher = sp.stats.fisher_exact([[A, B], [C, D]])
    print "##   P-value (Fisher):", p_value_fisher

    p_value_lrt, p_value_wald = binomial_logistic(A, B, C, D, num_parts)

    print "##   P-value (LRT):", p_value_lrt
    print "##   P-value (Wald):", p_value_wald

    return p_value_pearson, p_value_fisher, p_value_lrt, p_value_wald

def binomial_logistic(A_total, B_total, C_total, D_total, num_parts):
    """Run Binomial Logistic Regression.

    Note
    ----
    Divide the counts into 'num_parts'. This help prevent Perfect Seperation issue.
    Note that currently there is a random component to this split process (see
    split_num_into_parts() function).
    """

    print "##"
    print "## Binomial Logistic Regressionl: logit(pi) = B0 + B1*X1"

    num_positions = 2 * num_parts

    A_count_so_far = 0
    B_count_so_far = 0
    C_count_so_far = 0
    D_count_so_far = 0

    ## Split counts to 'num_parts' to revent prefect seperation.
    A_parts = split_num_into_parts(A_total, num_parts)
    B_parts = split_num_into_parts(B_total, num_parts)
    C_parts = split_num_into_parts(C_total, num_parts)
    D_parts = split_num_into_parts(D_total, num_parts)

    response_list = list()
    explanatory_list = list()

    assert sum(A_parts) == A_total
    assert sum(B_parts) == B_total
    assert sum(C_parts) == C_total
    assert sum(D_parts) == D_total

    for index in range(0, num_parts):
        A_count = A_parts[index]
        B_count = B_parts[index]

        explanatory_list.append([1.0, 1.0])
        response_list.append([A_count, B_count])

        A_count_so_far += A_count
        B_count_so_far += B_count

    for index in range(0, num_parts):
        C_count = C_parts[index]
        D_count = D_parts[index]

        explanatory_list.append([0.0, 1.0])
        response_list.append([C_count, D_count])

        C_count_so_far += C_count
        D_count_so_far += D_count

    # Consistency checks
    assert A_count_so_far == A_total
    assert B_count_so_far == B_total
    assert C_count_so_far == C_total
    assert D_count_so_far == D_total
    assert len(explanatory_list) == 2 * num_parts
    assert len(response_list) == 2 * num_parts

    # for index in range(len(response_list)):
    #     print "## %s --> %s" % (explanatory_list[index], response_list[index])

    glm_binom = sm.GLM(np.asarray(response_list), np.asarray(explanatory_list), family=sm.families.Binomial())
    fit_result = glm_binom.fit()

    print fit_result.summary()

    likelihood_ratio_test_statistic = 2 * (fit_result.llf - fit_result.llnull)
    # Survival function (also defined as 1 - cdf, but sf can handle very small p-values).
    p_value_lrt = sp.stats.chi2.sf(likelihood_ratio_test_statistic, 1)

    p_value_wald = fit_result.pvalues[0]

    print "##   llf:",  fit_result.llf
    print "##   llnull:",  fit_result.llnull
    print "##   likelihood_ratio_test_statistic:", likelihood_ratio_test_statistic

    return p_value_lrt, p_value_wald

def split_num_into_parts(total, k):
    """Split the total number into k parts (somewhat randomly)."""

    while True:
        sum_so_far = 0
        return_list = list()

        for index in range(0, k-1):

            avg_remaining_val = (total - sum_so_far) / (k - index)

            val = random.randint(0, (avg_remaining_val * 2))
            sum_so_far += val
            return_list.append(val)

        if sum_so_far <= total:
            return_list.append(total - sum_so_far)
            assert len(return_list) == k

            return return_list

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--dataset_name", action="store", required=True,
                        help="Either 'METABRIC' or 'TOY'.")

    parser.add_argument("--out_table", action="store", required=True,
                        metavar='FILE',
                        help="Path to output table")

    parser.add_argument("--out_plot", action="store", required=True,
                        metavar='FILE',
                        help="Path to output plot.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   dataset_name:", repr(options.dataset_name)
    print "##   out_table:", repr(options.out_table)
    print "##   out_plot", repr(options.out_plot)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
