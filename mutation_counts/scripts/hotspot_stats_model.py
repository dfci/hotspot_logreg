#!/usr/bin/env python2

"""
SYNOPSIS

    Statistically detects codon-level mutational hotspots using a Binomial
    Logistic Regression Model (with Fixed Effects).

WARNING

    This script contains PROTOTYPE CODE. Need to significantly test, refactor,
    clean-up and document this code before Publication.

EXAMPLES

    # (1) Analyze only DFCI and METABRIC dataset using Realistic Logistic Regression Model.
    ./hotspot_stats_model.py \
        --in_maf private/dfci_plus_metabric.r6.filtered_clinical.vep.maf \
        --in_map  private/extract_sequence/dfci_plus_metabric.r6.filtered_clinical.codon2trinuc_map.v1.txt \
        --in_source_maf_ids dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-2.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-3.r4.vep.maf \
                            metabric.vep.maf \
        --in_target_beds ../examples/inputs/private/dfci.PROFILE_POPv1.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv2.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv3.r1.bed \
                         ../examples/inputs/public/metabric.targetedIntervals.r1.bed \
        --out_table private/hotspot_stats_model_results/dfci_plus_metabric.hotspots.dev.1.txt \
        > private/hotspot_stats_model_results/log1.dfci_plus_metabric.dev.txt

    # (2) Analyze only DFCI and METABRIC dataset using a Simplified (Toy) Logistic
    #     Regression Model (i.e. No Fixed Effects).
    ./hotspot_stats_model.py \
        --in_maf private/dfci_plus_metabric.r6.filtered_clinical.vep.maf \
        --in_map  private/extract_sequence/dfci_plus_metabric.r6.filtered_clinical.codon2trinuc_map.v1.txt \
        --in_source_maf_ids dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-2.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-3.r4.vep.maf \
                            metabric.vep.maf \
        --in_target_beds ../examples/inputs/private/dfci.PROFILE_POPv1.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv2.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv3.r1.bed \
                         ../examples/inputs/public/metabric.targetedIntervals.r1.bed \
        --use_toy_stats_model \
        --out_table private/hotspot_stats_model_results/dfci_plus_metabric.hotspots.toy.1.txt \
        > private/hotspot_stats_model_results/log1.dfci_plus_metabric.toy.txt

    # (3) Analyze DFCI + Public dataset
    ./hotspot_stats_model.py \
        --in_maf private/merged.dfci_plus_public.r6.filtered_clinical.vep.maf \
        --in_map private/extract_sequence/merged.dfci_plus_public.r6.filtered_clinical.codon2trinuc_map.r1.txt \
        --in_source_maf_ids dfci.DFCI-ONCOPANEL-1.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-2.r4.vep.maf \
                            dfci.DFCI-ONCOPANEL-3.r4.vep.maf \
                            genie_msk.IMPACT341.r3.vep.maf \
                            genie_msk.IMPACT410.r3.vep.maf \
                            metabric.vep.maf \
                            sanger_somatic_coding_annotated.maf \
                            tcga.cell_2015.r2.vep.maf \
        --in_target_beds ../examples/inputs/private/dfci.PROFILE_POPv1.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv2.filtered.r1.bed \
                         ../examples/inputs/private/dfci.PROFILE_POPv3.r1.bed \
                         ../examples/inputs/public/genie/genie_msk.IMPACT341.r1.bed \
                         ../examples/inputs/public/genie/genie_msk.IMPACT410.r1.bed \
                         ../examples/inputs/public/metabric.targetedIntervals.r1.bed \
                         ../examples/inputs/public/from_intelccc/sanger_target_r1.bed \
                         ../examples/inputs/public/tcga_brca/tcga_target_v2.bed \
        --out_table private/hotspot_stats_model_results/merged.dfci_plus_public.hotspots.dev.11.txt \
        > private/hotspot_stats_model_results/log10.merged.dfci_plus_public.dev.txt

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import argparse
import math
import pandas as pd
import numpy as np
import scipy as sp
import traceback
import random
import itertools

import statsmodels.api as sm
import statsmodels.sandbox.stats.multicomp
import pybedtools

from collections import namedtuple

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 1000)

VERSION = 'mc_v0.9_dev'

LETTERS = ('A', 'C', 'G', 'T')

GenomicCoordID = namedtuple("GenomicCoordID", ["chrom", "position"])
CodonID = namedtuple("CodonID", ["gene", "codon_pos"])

INCLUDE_TRINUCLEOTIDE_FIXED_EFFECT = True
INCLUDE_GENE_FIXED_EFFECT = False
INCLUDE_HYPERMUTATED_FIXED_EFFECT = True

EXCLUDE_PUTATIVE_HOTSPOTS_FROM_BACKGROUND = True
EXCLUDE_FROM_BACKGROUND_COUNT_CUTOFF = 5

HIDE_SENSITIVE_LOGS = True

MIN_MUTATIONS_PER_GENE = 0

ANALYZE_ONLY_TEST_CODON = False

CENTERS_WITH_PANELS = [
    'DFCI',
    'OICR',
    'OHSU',
    'MSK-IMPACT341',
    'MSK-IMPACT410',
    'METABRIC',
    'GRCC',
    'VICC',
    'UHN',
    'MDA',
    'FOUNDATION'
]

def main(options):

    maf_df = import_annotated_maf(options.in_maf)

    # Note
    # ----
    # Important to tally-up samples before filtering MAF for Missense SNP variant
    # calls since samples with no Missense SNP variant calls will disappear after
    # the filtering.
    #
    # Import information about how many samples are targeted/not-targeted
    # at each BED interval.
    source_maf_id2target_intervals = map_source_maf_id_to_target_intervals(options.in_source_maf_ids,
                                                                           options.in_target_beds)

    source_maf_id2num_samples = map_source_maf_id_to_hardcoded_num_samples(options.in_source_maf_ids, maf_df)

    source_maf_id2center = map_source_maf_id_to_center(options.in_source_maf_ids, maf_df)

    maf_df = apply_missense_snp_filter(maf_df)

    maf_df = apply_exac_filter(maf_df)

    prefiltered_genes = get_prefiltered_genes(maf_df)

    maf_df = maf_df[maf_df['Gene_Symbol'].isin(prefiltered_genes)]

    transcriptcodon2trinucleotides = import_transcriptcodon2trinuc_map(options.in_map, prefiltered_genes)

    out_data_list = run_wrapper(maf_df,
                                transcriptcodon2trinucleotides,
                                source_maf_id2target_intervals,
                                source_maf_id2num_samples,
                                source_maf_id2center,
                                options.min_count,
                                VERSION,
                                options.use_toy_stats_model)

    if not ANALYZE_ONLY_TEST_CODON:
        write_output_table(out_data_list, options.out_table)


def run_wrapper(maf_df,
                transcriptcodon2trinucleotides,
                source_maf_id2target_intervals,
                source_maf_id2num_samples,
                source_maf_id2center,
                min_count,
                version_check,
                use_toy_stats_model=False):
    """Wrapper function to perform Statistical Analyzes detects codon-level mutational hotspots
    using a Binomial Logistic Regression Model (with Fixed Effects)."""

    if version_check != VERSION:
        # This check is to ensure that all deployed scripts inside IntelCCC;s Reliance Point are pulled from the same
        # version of the Git repo. Very unelegant solution but necessary due to limitation of Reliance Point.
        print "##"
        print "## ERROR: Detected Version Incompatibility inside %s's run_wrapper function." % os.path.basename(__file__)
        print "## ERROR: Specified version_check:", repr(version_check)
        print "## ERROR: Script's VERSION:", repr(VERSION)
        raise Exception("Detected Version Incompatibility inside %s's run_wrapper function." % os.path.basename(__file__))

    print "##", "-" * 50
    print "## Global Parameters:"
    print "##   INCLUDE_TRINUCLEOTIDE_FIXED_EFFECT:", repr(INCLUDE_TRINUCLEOTIDE_FIXED_EFFECT)
    print "##   INCLUDE_GENE_FIXED_EFFECT:", repr(INCLUDE_GENE_FIXED_EFFECT)
    print "##   INCLUDE_HYPERMUTATED_FIXED_EFFECT:", repr(INCLUDE_HYPERMUTATED_FIXED_EFFECT)
    print "##   EXCLUDE_PUTATIVE_HOTSPOTS_FROM_BACKGROUND:", repr(EXCLUDE_PUTATIVE_HOTSPOTS_FROM_BACKGROUND)
    print "##   EXCLUDE_FROM_BACKGROUND_COUNT_CUTOFF:", repr(EXCLUDE_FROM_BACKGROUND_COUNT_CUTOFF)
    print "##   MIN_MUTATIONS_PER_GENE:", repr(MIN_MUTATIONS_PER_GENE)
    print "##   ANALYZE_ONLY_TEST_CODON:", repr(ANALYZE_ONLY_TEST_CODON)
    print "##", "-" * 50

    print "##", "-" * 50
    print "## Function Parameters:"
    print "##   min_count:", repr(min_count)
    print "##   use_toy_stats_model:", repr(use_toy_stats_model)
    print "##", "-" * 50

    gene2transcript = get_gene2transcript(maf_df)

    supported_transcripts = set(transcriptcodon2trinucleotides.keys())

    num_samples = sum(source_maf_id2num_samples[key]['Total'] for key in source_maf_id2num_samples)

    genes = get_analyze_gene_list(maf_df, gene2transcript, supported_transcripts)

    print "##"
    print "## DEBUG: Num_Samples = %s" % num_samples
    print "## DEBUG: Num_Genes = %s" % len(genes)

    hypermutated_samples = get_hypermutated_samples(maf_df, source_maf_id2num_samples)

    genomiccoord2num_mutations = get_genomiccoord2num_mutations(maf_df, hypermutated_samples)

    genomiccoord2num_targeted = get_genomiccoord2num_targeted(genes,
                                                              gene2transcript,
                                                              transcriptcodon2trinucleotides,
                                                              source_maf_id2target_intervals,
                                                              source_maf_id2center,
                                                              source_maf_id2num_samples,
                                                              genomiccoord2num_mutations)

    codon_id2num_mutations = get_codon_id2num_mutations(genes,
                                                        gene2transcript,
                                                        transcriptcodon2trinucleotides,
                                                        genomiccoord2num_mutations)

    codon_id2exclude_from_background = get_codon_id2exclude_from_background(genes,
                                                                            gene2transcript,
                                                                            transcriptcodon2trinucleotides,
                                                                            codon_id2num_mutations)

    genes = sort_genes_by_mutability(genes,
                                     gene2transcript,
                                     transcriptcodon2trinucleotides,
                                     codon_id2num_mutations,
                                     codon_id2exclude_from_background)

    if not use_toy_stats_model:
        (fixedEffectVars2reponses,
         fixedEffectVarNames) = create_fixed_effect_vars2reponses_dict(genes,
                                                                       gene2transcript,
                                                                       transcriptcodon2trinucleotides,
                                                                       genomiccoord2num_mutations,
                                                                       codon_id2num_mutations,
                                                                       codon_id2exclude_from_background,
                                                                       genomiccoord2num_targeted)
    else:
        (fixedEffectVars2reponses,
         fixedEffectVarNames) = create_toy_fixed_effect_vars2reponses_dict(genes,
                                                                           gene2transcript,
                                                                           transcriptcodon2trinucleotides,
                                                                           genomiccoord2num_mutations,
                                                                           codon_id2num_mutations,
                                                                           codon_id2exclude_from_background,
                                                                           genomiccoord2num_targeted)

    # Analyze a specific test gene+codon of interest
    run_logistic_regression_specific_gene_codon(fixedEffectVars2reponses,
                                                fixedEffectVarNames,
                                                codon_id2num_mutations,
                                                "PIK3CA", 93)

    if ANALYZE_ONLY_TEST_CODON:
        # Finish analyzing the above test codon. Exit without analyzing the other
        # codons.
        return None

    out_data_list = run_logistic_regression_all_codons(fixedEffectVars2reponses,
                                                       fixedEffectVarNames,
                                                       codon_id2num_mutations,
                                                       gene2transcript,
                                                       min_count)

    return out_data_list

def get_prefiltered_genes(maf_df):
    """Perform prefiltering on the analyze gene list. Only include the gene
    if it appear in at least one of the sequencing panels (not Exome or WGS).

    TODO
    ----
    Expose command-line option to specify which source_maf_id correspond to a
    'sequencing panel' instead of hard-coding them inside the function.

    Alternatively, can also infer this from number of targeted genes. Perhaps
    define a panel as any source_maf_id that have lesser than 500 genes.
    """

    print "##"
    print "## Perform Prefiltering of Stat Model Analyze Gene List:"

    all_genes = sorted(maf_df['Gene_Symbol'].unique())
    gene2num_mutations = get_gene2num_mutations(maf_df)
    assert set(all_genes) == set(gene2num_mutations)

    # Get Gene that are part of a Panel (i.e. not Exome or WGS).
    center2genes = dict()
    for index, row in maf_df.iterrows():

        gene = row['Gene_Symbol']
        center = row['Center']

        if center not in center2genes:
            center2genes[center] = set()

        center2genes[center].add(gene)

    print "##"
    for center in center2genes:
        print "## DEBUG: '%s' dataset contains somatic mutation(s) at %s genes." % (center, len(center2genes[center]))

    panel_genes = set()
    for center in CENTERS_WITH_PANELS:
        if center in center2genes:
            panel_genes.update(center2genes[center])

    filtered_genes = panel_genes

    print "##"
    print "## DEBUG: Num_Genes (Before Prefiltering): %s" % len(all_genes)
    print "##"
    print "## DEBUG: Num_Genes (Part of a Panel): %s" % len(panel_genes)
    print "##"
    print "## DEBUG: Num_Genes (After Prefiltering): %s" % len(filtered_genes)

    return filtered_genes

def write_output_table(out_data_list, outfile):
    """Write the data in the out_data_list to the specified outfile."""

    out_header_cols = ['Gene_Symbol', 'Transcript_ID', 'Codon_Pos',
                       'Num_Samples_with_Mutation', 'PValue', 'QValue',
                       'LogLikelihoodRatio', 'LogLikelihood_Alt', 'LogLikelihood_Null']

    fout = open(outfile, 'w')
    fout.write("\t".join(out_header_cols) + '\n')

    for out_data in out_data_list:
        cols = list()
        cols.append(out_data['Gene_Symbol'])
        cols.append(out_data['Transcript_ID'])
        cols.append(str(out_data['Codon_Pos']))
        cols.append(str(out_data['Num_Samples_with_Mutation']))
        cols.append("%.3e" % out_data['PValue'])
        cols.append("%.3e" % out_data['QValue'])
        cols.append("%.3f" % out_data['LogLikelihoodRatio'])
        cols.append("%.3f" % out_data['LogLikelihood_Alt'])
        cols.append("%.3f" % out_data['LogLikelihood_Null'])
        fout.write("\t".join(cols) + '\n')

    fout.close()


def run_logistic_regression_all_codons(fixedEffectVars2reponses,
                                       fixedEffectVarNames,
                                       codon_id2num_mutations,
                                       gene2transcript,
                                       min_count):
    """
    Run Logistic Regression on all codon positions with mutation count equal
    or above the specified min_count threshold.
    """

    sorted_codon_ids = sorted(codon_id2num_mutations.iterkeys(), key=lambda(k): (k.gene, k.codon_pos))

    total_enumerated_codons = 0
    total_evaluated_codons = 0

    enumerated_genes = set()
    evaluated_genes = set()

    out_data_list = list()

    for codon_id in sorted_codon_ids:

        total_enumerated_codons += 1
        enumerated_genes.add(codon_id.gene)
        mutation_count = codon_id2num_mutations[codon_id]

        if mutation_count >= min_count:
            total_evaluated_codons += 1
            evaluated_genes.add(codon_id.gene)

            print "##"
            print "## [Putative Hotspot #%04d] Run Logistic Regression on Null Model [Gene: %s" % (total_evaluated_codons, codon_id.gene),
            if not HIDE_SENSITIVE_LOGS:
                print ",Codon_Pos=%s, Mutation_Count=%s" % (codon_id.codon_pos, mutation_count),
            print "]:"
            log_likelihood_null = run_binomail_regression(fixedEffectVars2reponses,
                                                          fixedEffectVarNames,
                                                          codon_id,
                                                          add_is_COI_indicator=False,
                                                          debug=False)

            print "##"
            print "## [Putative Hotspot #%04d] Run Logistic Regression on Alt Model [Gene: %s" % (total_evaluated_codons, codon_id.gene),
            if not HIDE_SENSITIVE_LOGS:
                print ",Codon_Pos=%s, Mutation_Count=%s" % (codon_id.codon_pos, mutation_count),
            print "]:"
            log_likelihood_alt = run_binomail_regression(fixedEffectVars2reponses,
                                                         fixedEffectVarNames,
                                                         codon_id,
                                                         add_is_COI_indicator=True,
                                                         debug=False)

            log_likelihood_ratio = log_likelihood_alt - log_likelihood_null
            llr_test_statistic = 2 * log_likelihood_ratio
            p_value_lrt = sp.stats.chi2.sf(llr_test_statistic, df=1)

            out_data = dict()
            out_data['Gene_Symbol'] = codon_id.gene
            out_data['Transcript_ID'] = gene2transcript[codon_id.gene]
            out_data['Codon_Pos'] = codon_id.codon_pos
            out_data['Num_Samples_with_Mutation'] = mutation_count
            out_data['PValue'] = p_value_lrt
            out_data['LogLikelihoodRatio'] = log_likelihood_ratio
            out_data['LogLikelihood_Alt'] = log_likelihood_alt
            out_data['LogLikelihood_Null'] = log_likelihood_null

            out_data_list.append(out_data)

        sys.stdout.flush()
        sys.stderr.flush()

    compute_q_values_wrapper(out_data_list, total_enumerated_codons)

    print "##"
    print "## ------------------------------------------------------------------"
    print "## Run Logistic Regression All Codons Summary:"
    print "##    Total Enumerated Codons: %s" % total_enumerated_codons
    print "##    Total Evaluated Codons: %s" % total_evaluated_codons
    print "##"
    print "##    Total Enumerated Genes: %s" % len(enumerated_genes)
    print "##    Total Evaluated Genes: %s" % len(evaluated_genes)
    print "## ------------------------------------------------------------------"

    return out_data_list

def run_logistic_regression_specific_gene_codon(fixedEffectVars2reponses,
                                                fixedEffectVarNames,
                                                codon_id2num_mutations,
                                                gene_of_interest,
                                                codon_of_interest):
    """Run Logistic Regression on a specific gene_codon

    Notes
    -----
    (1) Use this function to quickly test the Regression Logistic Code.

    (2) Here are some good test examples:

        # Low-count 'hotspot' in 'Hot' Gene (3 counts in DFCI + Public)
        (A) gene_of_interest = "PIK3CA" | codon_of_interest = 93

        # Medium-count 'hotspot' in 'Hot' Gene (18 counts in DFCI + Public)
        (B) gene_of_interest = "PIK3CA" | codon_of_interest = 118

        # Low-count 'hotspot' in 'Cold' Gene (3 counts in DFCI + Public)
        (C) gene_of_interest = "FANCF" | codon_of_interest = 216

        # Medium-count 'hotspot' in 'Hot' Gene (3 counts in DFCI + Public)
        (D) gene_of_interest = "TP53" | codon_of_interest = 130
    """

    codon_id_of_interest = CodonID(gene=gene_of_interest, codon_pos=codon_of_interest)

    print "##"
    print "## Run Logistic Regression on Null Model (Fixed Effect Vars Only):"
    log_likelihood_null = run_binomail_regression(fixedEffectVars2reponses,
                                                  fixedEffectVarNames,
                                                  codon_id_of_interest,
                                                  add_is_COI_indicator = False,
                                                  debug = True)

    print "##"
    print "## Run Logistic Regression on Alt Model [Add Codon of Interest Indicator]:"
    log_likelihood_alt = run_binomail_regression(fixedEffectVars2reponses,
                                                 fixedEffectVarNames,
                                                 codon_id_of_interest,
                                                 add_is_COI_indicator = True,
                                                 debug = True)

    log_likelihood_ratio = log_likelihood_alt - log_likelihood_null
    llr_test_statistic = 2 * log_likelihood_ratio

    p_value_lrt = sp.stats.chi2.sf(llr_test_statistic, df=1)

    print "##"
    print "## ------------------------------------------------------------------"
    print "## Logistic Regression Summary:"
    print "##    Codon of Interest: %s" % str(codon_id_of_interest)
    print "##    Mutation_Count: %s" % codon_id2num_mutations[codon_id_of_interest]
    print "##    log_likelihood_null: %s" % log_likelihood_null
    print "##    log_likelihood_alt: %s" % log_likelihood_alt
    print "##    log_likelihood_ratio: %s" % log_likelihood_ratio
    print "##    p_value_lrt: %s" % p_value_lrt
    print "## ------------------------------------------------------------------"

def run_binomail_regression(fixedEffectVars2reponses,
                            fixedEffectVarNames,
                            codon_id_of_interest,
                            add_is_COI_indicator = True,
                            debug = False):
    """Run the Binomail_Regression Model and return the Log-Likelihood value."""


    assert codon_id_of_interest != None

    start_time = time.time()

    attempt_num = 0

    log_likelihood = None  # Return Variable

    while True:
        attempt_num += 1
        num_parts = attempt_num

        try:
            # Construct explanatory and response matrix
            explanatory_list, response_list = construct_explanatory_and_response_matrix(fixedEffectVars2reponses,
                                                                                        codon_id_of_interest,
                                                                                        add_is_COI_indicator,
                                                                                        num_parts)

            start_time2 = time.time()
            glm_binom = sm.GLM(np.asarray(response_list),
                               np.asarray(explanatory_list),
                               family=sm.families.Binomial())
            fit_result = glm_binom.fit()
            log_likelihood = fit_result.llf - fit_result.llnull

            if debug:
                explanatory_var_names = list()
                print "##"
                print "## Time taken (Fit Binomail Regression): %.3f secs" % (time.time() - start_time2)
                if add_is_COI_indicator:
                    is_gene_codon_of_interest = 'Is_%s_%s' % (codon_id_of_interest.gene, codon_id_of_interest.codon_pos)
                    explanatory_var_names = [is_gene_codon_of_interest] + fixedEffectVarNames
                else:
                    explanatory_var_names = fixedEffectVarNames
                print "##"
                print "## Logistic Regression Summary:"
                print fit_result.summary(xname=explanatory_var_names)

        except:
            print "##"
            print "## WARNING: Logistic Regression failed to converge for attempt #%d (num_parts = %d)." % (attempt_num, num_parts),
            print "Here is the error message:"
            print "##"
            traceback.print_exc()
            print "##"
            print "## WARNING: Retrying ..."
            if attempt_num == 1:
                raise Exception("Logistic Regression failed to converge after %d attempts!" % attempt_num)
        else:
            if attempt_num > 1:
                print "##"
                print "## Successfully Ran Logistic Regression for after multiple attempts (#%d) (num_parts = %d)." % (attempt_num, num_parts)
            break

    print "##"
    print "## Time taken (Construct Array + Run Binomail Regression): %.3f secs" % (time.time() - start_time)

    return log_likelihood

def construct_explanatory_and_response_matrix(fixedEffectVars2reponses,
                                              codon_id_of_interest,
                                              add_is_COI_indicator,
                                              num_parts):
    """Create Explanatory and Response Matrices to be input to the Binomail
    Logistic Regression Model.

    Note
    ----
    (1) If add_is_COI_indicator is True, then include 'is_COI' indicator as one
        of the explanatory variable (which indicate if the data correspond to the
        codon-of-interest-position.


    (2) For each fixed_effect_vars, divide the response counts into 'num_parts'.
        This help prevent Perfect Seperation issue. Note that currently there is
        a random component to this division (see split_num_into_parts() function).
    """

    assert codon_id_of_interest != None

    explanatory_list = list()
    response_list = list()

    for fixed_effect_vars, response_vars in fixedEffectVars2reponses.iteritems():

        total_mut_count = response_vars['Total']['Mut_Count']
        total_ref_count = response_vars['Total']['Ref_Count']

        # Divide into parts to help prevent Perfect Seperation issue.
        total_mut_count_parts = split_num_into_parts(total_mut_count, num_parts)
        total_ref_count_parts = split_num_into_parts(total_ref_count, num_parts)

        if codon_id_of_interest in response_vars:

            # coi: Codon-Of-Interest
            ref_count_at_coi = response_vars[codon_id_of_interest]['Ref_Count']
            mut_count_at_coi = response_vars[codon_id_of_interest]['Mut_Count']
            added_cio_to_total = response_vars[codon_id_of_interest]['Added_To_Total']

            ref_count_at_coi_parts = split_num_into_parts(ref_count_at_coi, num_parts)
            mut_count_at_coi_parts = split_num_into_parts(mut_count_at_coi, num_parts)

            if added_cio_to_total:
                # Codon-Of-Interest was not excluded from total_counts, so
                # need adjustment.
                ref_count_at_other_parts = split_num_into_parts(total_ref_count - ref_count_at_coi, num_parts)
                mut_count_at_other_parts = split_num_into_parts(total_mut_count - mut_count_at_coi, num_parts)
            else:
                # Codon-Of-Interest was already excluded from total_counts,
                # so no adjustment necessary.
                ref_count_at_other_parts = split_num_into_parts(total_ref_count, num_parts)
                mut_count_at_other_parts = split_num_into_parts(total_mut_count, num_parts)

            for index in range(num_parts):
                if add_is_COI_indicator:
                    # Include explanatory variable for the codon_id_of_interest in the
                    # Regression Model.
                    explanatory_list.append([0] + list(fixed_effect_vars))
                    explanatory_list.append([1] + list(fixed_effect_vars))
                else:
                    explanatory_list.append(list(fixed_effect_vars))
                    explanatory_list.append(list(fixed_effect_vars))

                response_list.append([mut_count_at_other_parts[index], ref_count_at_other_parts[index]])
                response_list.append([mut_count_at_coi_parts[index], ref_count_at_coi_parts[index]])

        else:
            for index in range(num_parts):
                if add_is_COI_indicator:
                    explanatory_list.append([0] + list(fixed_effect_vars))
                else:
                    explanatory_list.append(list(fixed_effect_vars))

                response_list.append([total_mut_count_parts[index], total_ref_count_parts[index]])

    assert len(explanatory_list) == len(response_list)

    return explanatory_list, response_list

def get_codon_id2exclude_from_background(genes,
                                         gene2transcript,
                                         transcriptcodon2trinucleotides,
                                         codon_id2num_mutations):
    """Create and return a dictionary which indicate whether the data at each
    codon_id should be excluded from background mutation rate computation.

    Notes
    -----
    (1) The idea is the background mutation rate may be over-estimated if we
        include hotspot positions in the ca. So we want 'rough' filter that will
        identify and exclude likely hotspots.
    """

    codon_id2exclude_from_background = dict()

    for gene in genes:
        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):
            codon_id = CodonID(gene=gene, codon_pos=codon_position)

            assert codon_id not in codon_id2exclude_from_background
            codon_id2exclude_from_background[codon_id] = False

            if EXCLUDE_PUTATIVE_HOTSPOTS_FROM_BACKGROUND:
                if codon_id2num_mutations[codon_id] >= EXCLUDE_FROM_BACKGROUND_COUNT_CUTOFF:
                    codon_id2exclude_from_background[codon_id] = True

    return codon_id2exclude_from_background

def get_codon_id2num_mutations(genes,
                               gene2transcript,
                               transcriptcodon2trinucleotides,
                               genomiccoord2num_mutations):
    """Create and return a dictionary mapping each codon_id to the number of
    mutation observed at the codon_id."""

    codon_id2num_mutations = dict()

    for gene in genes:
        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):
            codon_id = CodonID(gene=gene, codon_pos=codon_position)

            assert codon_id not in codon_id2num_mutations
            codon_id2num_mutations[codon_id] = 0

            for frame in [0,1,2]:
                data = transcriptcodon2trinucleotides[transcript][codon_position][frame]
                genomic_coord_id = GenomicCoordID(chrom=data['Chrom'], position=int(data['Position']))

                if genomic_coord_id in genomiccoord2num_mutations:
                    num_sample_with_mutation = genomiccoord2num_mutations[genomic_coord_id]['Total']
                else:
                    num_sample_with_mutation = 0

                codon_id2num_mutations[codon_id] += num_sample_with_mutation

    return codon_id2num_mutations

def create_toy_fixed_effect_vars2reponses_dict(genes,
                                               gene2transcript,
                                               transcriptcodon2trinucleotides,
                                               genomiccoord2num_mutations,
                                               codon_id2num_mutations,
                                               codon_id2exclude_from_background,
                                               genomiccoord2num_targeted):
    """Create a dictionary mapping each fixed_effect_vars combination to
    observed mutation counts and observed ref counts.

    Notes
    -----
    This is not Toy (Simplified) Version with only 'Intercept' (no actual
    fixed effects).
    """

    print "##"
    print "## Enter create_toy_fixed_effect_vars2reponses_dict()"

    start_time = time.time()

    fixedEffectVars2reponses = dict()
    fixedEffectVarNames = ['Intercept']

    total_mutation_events = 0
    for gene_index, gene in enumerate(genes, 1):

        print "##"
        print "## PROGRESS: [%03d] %-9s --> create_fixed_effect_vars2reponses_dict" % (gene_index, gene)

        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):

            for frame in [0,1,2]:

                data = transcriptcodon2trinucleotides[transcript][codon_position][frame]
                genomic_coord_id = GenomicCoordID(chrom=data['Chrom'], position=int(data['Position']))
                codon_id = CodonID(gene=gene, codon_pos=codon_position)

                num_sample_targeted = genomiccoord2num_targeted[genomic_coord_id]['Total']

                # Get Fixed Effect Explanatory Tuple
                fixed_effect_vars = list()

                # Beta0 (Intercept) variable
                fixed_effect_vars.append(1)

                fixed_effect_vars = tuple(fixed_effect_vars)

                assert len(fixed_effect_vars) == len(fixedEffectVarNames)

                if fixed_effect_vars not in fixedEffectVars2reponses:
                    fixedEffectVars2reponses[fixed_effect_vars] = dict()
                    fixedEffectVars2reponses[fixed_effect_vars]['Total'] = dict()
                    fixedEffectVars2reponses[fixed_effect_vars]['Total']['Mut_Count'] = 0
                    fixedEffectVars2reponses[fixed_effect_vars]['Total']['Ref_Count'] = 0

                if codon_id not in fixedEffectVars2reponses[fixed_effect_vars]:
                    fixedEffectVars2reponses[fixed_effect_vars][codon_id] = dict()
                    fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Mut_Count'] = 0
                    fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Ref_Count'] = 0

                # Get Response Variable Counts.
                if genomic_coord_id in genomiccoord2num_mutations:
                    num_sample_with_mutation = genomiccoord2num_mutations[genomic_coord_id]['Total']
                else:
                    num_sample_with_mutation = 0

                assert num_sample_targeted >= num_sample_with_mutation

                num_sample_with_ref = num_sample_targeted - num_sample_with_mutation

                added_to_total = True
                if codon_id2exclude_from_background[codon_id]:
                    added_to_total = False

                if added_to_total:
                    fixedEffectVars2reponses[fixed_effect_vars]['Total']['Mut_Count'] += num_sample_with_mutation
                    fixedEffectVars2reponses[fixed_effect_vars]['Total']['Ref_Count'] += num_sample_with_ref

                if 'Added_To_Total' not in fixedEffectVars2reponses[fixed_effect_vars][codon_id]:
                    fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Added_To_Total'] = added_to_total
                else:
                    # Consistency Check
                    assert fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Added_To_Total'] == added_to_total

                fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Mut_Count'] += num_sample_with_mutation
                fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Ref_Count'] += num_sample_with_ref

                # DEBUG
                total_mutation_events += num_sample_with_mutation

    print "##"
    print "##   Total Mutation Events = %s." % total_mutation_events
    print "##"
    print "##   Time taken (Construct fixedEffectVars2reponses dict): %.3f secs" % (time.time() - start_time)

    return fixedEffectVars2reponses, fixedEffectVarNames

def create_fixed_effect_vars2reponses_dict(genes,
                                           gene2transcript,
                                           transcriptcodon2trinucleotides,
                                           genomiccoord2num_mutations,
                                           codon_id2num_mutations,
                                           codon_id2exclude_from_background,
                                           genomiccoord2num_targeted):
    """Create a dictionary mapping each fixed_effect_vars combination to
    observed mutation counts and observed ref counts."""

    print "##"
    print "## Enter create_fixed_effect_vars2reponses_dict()"

    start_time = time.time()

    ref_trinucleotides = get_ref_trinucleotides(genes,
                                                gene2transcript,
                                                transcriptcodon2trinucleotides,
                                                genomiccoord2num_mutations,
                                                codon_id2num_mutations,
                                                codon_id2exclude_from_background,
                                                genomiccoord2num_targeted)

    fixedEffectVars2reponses = dict()
    fixedEffectVarNames = ['Intercept']

    if INCLUDE_TRINUCLEOTIDE_FIXED_EFFECT:
        for index, key_id in enumerate(ref_trinucleotides):
            if index == 0:
                # Omit the first 'ref_trinucleotide' indicator term to avoid
                # multicollinearity with the constant (intercept) term.
                continue

            fixedEffectVarNames.append("Is_%s" % (key_id.ref_trinucleotide))

    if INCLUDE_GENE_FIXED_EFFECT:
        for fe_gene_index, fe_gene in enumerate(genes):
            if fe_gene_index == 0:
                # Omit the first 'gene' indicator term to avoid multicollinearity
                # with the constant (intercept) term.
                continue

            fixedEffectVarNames.append('Is_%s' % fe_gene)

    if INCLUDE_HYPERMUTATED_FIXED_EFFECT:
        fixedEffectVarNames.append('Is_Hypermutated')

    total_mutation_events = 0
    for gene_index, gene in enumerate(genes, 1):

        print "##"
        print "## PROGRESS: [%03d] %-9s --> create_fixed_effect_vars2reponses_dict" % (gene_index, gene)

        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):
            for frame in [0,1,2]:

                data = transcriptcodon2trinucleotides[transcript][codon_position][frame]
                genomic_coord_id = GenomicCoordID(chrom=data['Chrom'], position=int(data['Position']))
                codon_id = CodonID(gene=gene, codon_pos=codon_position)

                ref_trinucleotide = data['Trinucleotide']
                ref_allele = ref_trinucleotide[1]  # the middle nucleotide

                for hypermutated_key in ['Hypermutated', 'NotHypermutated']:

                    num_sample_targeted = genomiccoord2num_targeted[genomic_coord_id][hypermutated_key]

                    # Get Fixed Effect Explanatory Tuple
                    fixed_effect_vars = list()

                    # Beta0 (Intercept) variable
                    fixed_effect_vars.append(1)

                    if INCLUDE_TRINUCLEOTIDE_FIXED_EFFECT:
                        for fe_index, fe_key_id in enumerate(ref_trinucleotides):
                            if fe_index == 0:
                                # Omit the first 'ref_trinucleotide' indicator term to avoid
                                # multicollinearity with the constant (intercept) term.
                                continue

                            if ref_trinucleotide == fe_key_id.ref_trinucleotide:
                                fixed_effect_vars.append(1)
                            else:
                                fixed_effect_vars.append(0)

                    if INCLUDE_GENE_FIXED_EFFECT:
                        for fe_gene_index, fe_gene in enumerate(genes):
                            if fe_gene_index == 0:
                                # Omit the first 'gene' indicator term to avoid multicollinearity
                                # with the constant (intercept) term.
                                continue

                            if gene == fe_gene:
                                fixed_effect_vars.append(1)
                            else:
                                fixed_effect_vars.append(0)

                    if INCLUDE_HYPERMUTATED_FIXED_EFFECT:
                        if hypermutated_key == 'Hypermutated':
                            fixed_effect_vars.append(1)
                        else:
                            fixed_effect_vars.append(0)

                    fixed_effect_vars = tuple(fixed_effect_vars)

                    assert len(fixed_effect_vars) == len(fixedEffectVarNames)

                    if fixed_effect_vars not in fixedEffectVars2reponses:
                        fixedEffectVars2reponses[fixed_effect_vars] = dict()
                        fixedEffectVars2reponses[fixed_effect_vars]['Total'] = dict()
                        fixedEffectVars2reponses[fixed_effect_vars]['Total']['Mut_Count'] = 0
                        fixedEffectVars2reponses[fixed_effect_vars]['Total']['Ref_Count'] = 0

                    if codon_id not in fixedEffectVars2reponses[fixed_effect_vars]:
                        fixedEffectVars2reponses[fixed_effect_vars][codon_id] = dict()
                        fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Mut_Count'] = 0
                        fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Ref_Count'] = 0

                    # Get Response Variable Counts.
                    if genomic_coord_id in genomiccoord2num_mutations:
                        num_mutations_dict = genomiccoord2num_mutations[genomic_coord_id]
                        genomic_ref_allele = num_mutations_dict['RefAllele']
                        num_sample_with_mutation = num_mutations_dict[hypermutated_key]

                        # Consistency check
                        if genomic_ref_allele not in [ref_allele, complement_nt_letter(ref_allele)]:
                            print "## ERROR: Problem genomic_coord_id: %s" % str(genomic_coord_id)
                            print "## ERROR: Problem ref_allele: %s" % ref_allele
                            print "## ERROR: Problem genomic_ref_allele: %s" % genomic_ref_allele
                            raise Exception("Inconsistency between ref_allele and genomic_ref_allele!")
                    else:
                        num_sample_with_mutation = 0

                    assert num_sample_targeted >= num_sample_with_mutation

                    num_sample_with_ref = num_sample_targeted - num_sample_with_mutation

                    added_to_total = True
                    if codon_id2exclude_from_background[codon_id]:
                        added_to_total = False

                    if added_to_total:
                        fixedEffectVars2reponses[fixed_effect_vars]['Total']['Mut_Count'] += num_sample_with_mutation
                        fixedEffectVars2reponses[fixed_effect_vars]['Total']['Ref_Count'] += num_sample_with_ref

                    if 'Added_To_Total' not in fixedEffectVars2reponses[fixed_effect_vars][codon_id]:
                        fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Added_To_Total'] = added_to_total
                    else:
                        # Consistency Check
                        assert fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Added_To_Total'] == added_to_total

                    fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Mut_Count'] += num_sample_with_mutation
                    fixedEffectVars2reponses[fixed_effect_vars][codon_id]['Ref_Count'] += num_sample_with_ref

                    # DEBUG
                    total_mutation_events += num_sample_with_mutation

    print "##"
    print "##   Total Mutation Events = %s." % total_mutation_events
    print "##"
    print "##   Time taken (Construct fixedEffectVars2reponses dict): %.3f secs" % (time.time() - start_time)

    return fixedEffectVars2reponses, fixedEffectVarNames

def get_ref_trinucleotides(genes,
                           gene2transcript,
                           transcriptcodon2trinucleotides,
                           genomiccoord2num_mutations,
                           codon_id2num_mutations,
                           codon_id2exclude_from_background,
                           genomiccoord2num_targeted):
    """Create a list of all 32 possible reference trinucletotide contexts.

    Notes
    -----
    (1)  Sort the trinucleotide contexts by the mutation rate (ratio) at the
         central nucleotide.
    """

    KeyID = namedtuple("KeyID", ["ref_trinucleotide"])
    key2mutcount = dict()  # Y=1 counts
    key2allcount = dict()  # Y=0 or Y=1 counts
    key2ratio = dict()
    key_ids = list()

    trinucleotides = set()
    for values in itertools.product(LETTERS, LETTERS, LETTERS):
        values = [str(value) for value in values]
        trinucleotide = "".join(values)
        trinucleotide = normalize_trinucleotide(trinucleotide)
        trinucleotides.add(trinucleotide)
    trinucleotides = sorted(trinucleotides)

    # Initialize the count + ratio dictionaries
    for trinucleotide in trinucleotides:
        ref_allele = trinucleotide[1]  # the middle nucleotide
        key_id = KeyID(ref_trinucleotide=trinucleotide)
        key2mutcount[key_id] = 0
        key2allcount[key_id] = 0
        key2ratio[key_id] = 0.0
        key_ids.append(key_id)

    for gene_index, gene in enumerate(genes, 1):
        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):

            codon_id = CodonID(gene=gene, codon_pos=codon_position)

            if codon_id2exclude_from_background[codon_id]:
                continue

            for frame in [0,1,2]:
                data = transcriptcodon2trinucleotides[transcript][codon_position][frame]
                genomic_coord_id = GenomicCoordID(chrom=data['Chrom'], position=int(data['Position']))
                ref_trinucleotide = data['Trinucleotide']
                assert len(ref_trinucleotide) == 3

                num_sample_targeted = genomiccoord2num_targeted[genomic_coord_id]['Total']

                # Get mutation and target counts with each trinucleotide context
                if genomic_coord_id in genomiccoord2num_mutations:
                    num_sample_with_mutation = genomiccoord2num_mutations[genomic_coord_id]['Total']

                else:
                    num_sample_with_mutation = 0

                key_id = KeyID(ref_trinucleotide=ref_trinucleotide)
                key2mutcount[key_id] += num_sample_with_mutation
                key2allcount[key_id] += num_sample_targeted

    for key_id in key_ids:
        key2ratio[key_id] = key2mutcount[key_id] / float(key2allcount[key_id])

    sorted_key_ids = sorted(key_ids, key=key2ratio.get)

    debug = True

    if debug:
        total_mut_count = 0
        total_all_count = 0
        print "##"
        print "## DEBUG: Mutation and Targeted Counts at each Trinucleotide Context:"
        for index, key_id in enumerate(sorted_key_ids, start=1):
            total_mut_count += key2mutcount[key_id]
            total_all_count += key2allcount[key_id]
            print "##  %2d  %s --> %5d     %10d     %.2E" % (index,
                                                             str(key_id),
                                                             key2mutcount[key_id],
                                                             key2allcount[key_id],
                                                             key2ratio[key_id])

        print "##"
        print "## DEBUG: Aggregated Mutation and Targeted Counts across all Trinucleotide Context:"
        print "##  All --> %5d     %10d     %.2E" % (total_mut_count,
                                                     total_all_count,
                                                     total_mut_count/float(total_all_count))


    # Consistency check.
    assert len(sorted_key_ids) == 32

    return sorted_key_ids

def get_analyze_gene_list(maf_df, gene2transcript, supported_transcripts):
    """Get list of genes to perform analysis on.

    Note
    ----
    (1) Exclude a gene from analysis if there are fewer than 'MIN_MUTATIONS_PER_GENE'
        mutation events observed in the entire gene across all samples

    (2) Exclude a gene from analysis if it doesn't map to any
        supported_transcripts
    """

    all_genes = sorted(maf_df['Gene_Symbol'].unique())

    gene2num_mutations = get_gene2num_mutations(maf_df)

    assert set(all_genes) == set(gene2transcript)
    assert set(all_genes) == set(gene2num_mutations)

    to_analyze_genes = list()

    num_too_few_mutations_genes = 0
    num_unsupported_transcript_genes = 0

    for gene in all_genes:

        if gene2num_mutations[gene] < MIN_MUTATIONS_PER_GENE:
            num_too_few_mutations_genes += 1
            continue

        transcript = gene2transcript[gene]

        if transcript not in supported_transcripts:
            print "##"
            print "## WARNING: Excluding gene '%s' since it maps to a unsupported" % gene,
            print "transcript '%s'." % transcript
            num_unsupported_transcript_genes += 1
            continue

        to_analyze_genes.append(gene)

    print "##"
    print "## WARNING: Excluded %d Gene(s) (Reason: Gene has fewer than" % num_too_few_mutations_genes,
    print "%s Missense Mutations called across all Samples + all Codon Positions)." % MIN_MUTATIONS_PER_GENE

    print "##"
    print "## WARNING: Excluded %d Gene(s) (Reason: Gene maps to unsupported Transcript)." % num_unsupported_transcript_genes

    return to_analyze_genes

def get_gene2num_mutations(maf_df):

    gene2num_mutations = dict()

    for index, row in maf_df.iterrows():

        gene = row['Gene_Symbol']
        if gene not in gene2num_mutations:
            gene2num_mutations[gene] = 0

        gene2num_mutations[gene] += 1

    return gene2num_mutations

def sort_genes_by_mutability(genes,
                             gene2transcript,
                             transcriptcodon2trinucleotides,
                             codon_id2num_mutations,
                             codon_id2exclude_from_background):

    gene2num_background_mutations = dict()

    for codon_id in codon_id2num_mutations:

        gene = codon_id.gene

        added_to_background = True

        if codon_id2exclude_from_background[codon_id]:
            added_to_background = False

        if gene not in gene2num_background_mutations:
            gene2num_background_mutations[gene] = 0

        if added_to_background:
            gene2num_background_mutations[gene] += codon_id2num_mutations[codon_id]

    assert set(genes) == set(gene2num_background_mutations)

    gene2mutability = dict()

    for gene in genes:
        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])
        num_background_mutations = gene2num_background_mutations[gene]

        gene2mutability[gene] = dict()
        gene2mutability[gene]['protein_length'] = protein_length
        gene2mutability[gene]['mutability'] = num_background_mutations/float(protein_length)

    sorted_genes = sorted(genes, key=lambda(gene):gene2mutability[gene]['mutability'])

    debug = True

    if debug:
        print "##"
        print "## List of '%s' Genes included in Analysis:" % len(genes)
        print "##"
        for index, gene in enumerate(sorted_genes, 1):
            print "##      [%03d] %-9s --> (Num_Background_Mutations: %s, Protein_Length: %s)" % (index, gene,
                                                                                                  gene2num_background_mutations[gene],
                                                                                                  gene2mutability[gene]['protein_length'])

    assert set(genes) == set(sorted_genes)

    return sorted_genes

def get_genomiccoord2num_mutations(maf_df, hypermutated_samples):
    """Figure out information about the num_samples with mutations at each genomics
    coordinate."""

    genomiccoord2num_mutations = dict()

    all_centers = list(maf_df['Center'].unique())
    all_source_maf_ids = list(maf_df['Source_Maf'].unique())
    GenomicCoordID = namedtuple("GenomicCoordID", ["chrom", "position"])

    for index, row in maf_df.iterrows():

        chrom = row['Chromosome']
        position = int(row['Start_Position'])
        ref_allele = row['Reference_Allele']
        mut_allele = row['Tumor_Seq_Allele2']
        sample_id = row['Tumor_Sample_Barcode']

        center = row['Center']
        source_maf_id = row['Source_Maf']

        assert row['Start_Position'] == row['End_Position']
        assert ref_allele in LETTERS
        assert mut_allele in LETTERS
        assert ref_allele != mut_allele

        genomic_coord_id = GenomicCoordID(chrom=chrom, position=position)

        if genomic_coord_id not in genomiccoord2num_mutations:
            genomiccoord2num_mutations[genomic_coord_id] = initialize_num_mutations_dict(ref_allele, all_centers, all_source_maf_ids)
        else:
            assert genomiccoord2num_mutations[genomic_coord_id]['RefAllele'] == ref_allele

        if sample_id in hypermutated_samples:
            hypermutated_key = 'Hypermutated'
        else:
            hypermutated_key = 'NotHypermutated'

        for key in ['Total', hypermutated_key]:
            genomiccoord2num_mutations[genomic_coord_id][key] += 1
            genomiccoord2num_mutations[genomic_coord_id]['PerCenter'][center][key]  += 1
            genomiccoord2num_mutations[genomic_coord_id]['PerSourceMafID'][source_maf_id][key]  += 1

    return genomiccoord2num_mutations

def get_genomiccoord2num_targeted(genes,
                                  gene2transcript,
                                  transcriptcodon2trinucleotides,
                                  source_maf_id2target_intervals,
                                  source_maf_id2center,
                                  source_maf_id2num_samples,
                                  genomiccoord2num_mutations):
    """Figure out information about the num_samples targeted at each genomics
    coordinate."""

    genomiccoord2num_targeted = dict()

    all_source_maf_ids = sorted(source_maf_id2target_intervals.keys())
    all_centers = sorted(set(source_maf_id2center.values()))

    for gene_index, gene in enumerate(genes, 1):
        transcript = gene2transcript[gene]
        protein_length = max(transcriptcodon2trinucleotides[transcript])

        for codon_position in xrange(1, protein_length+1):
            for frame in [0,1,2]:
                data = transcriptcodon2trinucleotides[transcript][codon_position][frame]
                genomic_coord_id = GenomicCoordID(chrom=data['Chrom'], position=int(data['Position']))

                num_targeted_dict = initialize_num_targeted_dict(all_centers, all_source_maf_ids)

                for source_maf_id in all_source_maf_ids:
                    target_intervals = source_maf_id2target_intervals[source_maf_id]
                    point_interval = create_point_interval(genomic_coord_id.chrom, genomic_coord_id.position)
                    center = source_maf_id2center[source_maf_id]

                    is_targeted_pos = target_intervals.any_hits(point_interval)

                    is_mutation_called_pos = False
                    if genomic_coord_id in genomiccoord2num_mutations:
                        if source_maf_id in genomiccoord2num_mutations[genomic_coord_id]['PerSourceMafID']:
                            if genomiccoord2num_mutations[genomic_coord_id]['PerSourceMafID'][source_maf_id]['Total'] > 0:
                                num_mutations_dict = genomiccoord2num_mutations[genomic_coord_id]
                                is_mutation_called_pos = True

                    if is_mutation_called_pos and not is_targeted_pos:
                        # Gracefully handle potential inconsistency between Mutation Called positions and Targeted regions.
                        is_targeted_pos = True
                        print ""
                        print "## WARNING: Found Mutation(s) called at an Untargeted Codon. Changing Target Status of this Codon to 'Targeted'!"
                        print "## WARNING:     Problem Gene_Codon: '%s_%s'" % (gene, codon_position)
                        print "## WARNING:     Problem Genomic_Coordinate: '%s:%s'" % (genomic_coord_id.chrom, genomic_coord_id.position)
                        print "## WARNING:     Problem Center: '%s'" % center
                        print "## WARNING:     Problem SourceMafID: '%s'" % source_maf_id
                        print "## WARNING:     Num_Samples_with_Mutation:"
                        print "## WARNING:         Total         : %s" % num_mutations_dict['Total']
                        print "## WARNING:         At_Center     : %s" % num_mutations_dict['PerCenter'][center]['Total']
                        print "## WARNING:         At_SourceMafID: %s" % num_mutations_dict['PerSourceMafID'][source_maf_id]['Total']

                    if is_targeted_pos:
                        for key in ['Total', 'Hypermutated', 'NotHypermutated']:
                            count = source_maf_id2num_samples[source_maf_id][key]
                            num_targeted_dict[key] += count
                            num_targeted_dict['PerCenter'][center][key] += count
                            num_targeted_dict['PerSourceMafID'][source_maf_id][key] += count

                # Consistency check and update return dict
                if genomic_coord_id not in genomiccoord2num_targeted:
                    genomiccoord2num_targeted[genomic_coord_id] = num_targeted_dict
                else:
                    # This else clause is reached with in the rare situation where the coding
                    # regions of 2 different genes overlap. An example of this is cis-natural
                    # antisense transcripts such as 'PAXIP1' and PAXP1-AS2. For more details,
                    # see: https://en.wikipedia.org/wiki/Cis-natural_antisense_transcript.
                    print "##"
                    print "## WARNING: Current genomic_coord_id already exist in genomiccoord2num_targeted."
                    print "## WARNING:     Problem (Gene, Codon_Pos, Frame): (%s, %s, %s)" % (gene, codon_position, frame)
                    print "## WARNING:     Problem Genomic_Coordinate: '%s:%s'" % (genomic_coord_id.chrom, genomic_coord_id.position)

    return genomiccoord2num_targeted

def initialize_num_mutations_dict(ref_allele, all_centers, all_source_maf_ids):

    # Check for duplicate values
    assert len(set(all_centers)) == len(all_centers)
    assert len(set(all_source_maf_ids)) == len(all_source_maf_ids)

    keys = ['Hypermutated', 'NotHypermutated', 'Total']

    init_dict = dict()
    init_dict = {key: 0 for key in keys}
    init_dict['RefAllele'] = ref_allele

    # PerCenter Information
    init_dict['PerCenter'] = dict()
    for center in all_centers:
        init_dict['PerCenter'][center] = {key: 0 for key in keys}

    # PerSourceMafID Information
    init_dict['PerSourceMafID'] = dict()
    for source_maf_id in all_source_maf_ids:
        init_dict['PerSourceMafID'][source_maf_id] = {key: 0 for key in keys}

    return init_dict

def initialize_num_targeted_dict(all_centers, all_source_maf_ids):

    # Check for duplicate values
    assert len(set(all_centers)) == len(all_centers)
    assert len(set(all_source_maf_ids)) == len(all_source_maf_ids)

    keys = ['Hypermutated', 'NotHypermutated', 'Total']

    init_dict = dict()
    init_dict = {key: 0 for key in keys}

    # PerCenter Information
    init_dict['PerCenter'] = dict()
    for center in all_centers:
        init_dict['PerCenter'][center] = {key: 0 for key in keys}

    # PerSourceMafID Information
    init_dict['PerSourceMafID'] = dict()
    for source_maf_id in all_source_maf_ids:
        init_dict['PerSourceMafID'][source_maf_id] = {key: 0 for key in keys}

    return init_dict

def get_gene2transcript(maf_df):
    """Create and return a mapping from gene to transcript. If a gene maps to
    multiple transcript, then select the transcripts with the most data (e.g.
    mutations).
    """

    gene2transcript2count = dict()

    for index, row in maf_df.iterrows():

        gene = row['Gene_Symbol']
        transcript = row['Transcript_ID']

        if gene not in gene2transcript2count:
            gene2transcript2count[gene] = dict()

        if transcript not in gene2transcript2count[gene]:
            gene2transcript2count[gene][transcript] = 0

        gene2transcript2count[gene][transcript] += 1

    gene2transcript = dict()

    for gene in gene2transcript2count:

        transcripts = gene2transcript2count[gene]

        best_transcript = None
        best_count = None

        for transcript in sorted(transcripts):

            count = gene2transcript2count[gene][transcript]

            if best_transcript == None or count > best_count:
                best_transcript = transcript
                best_count = count

        gene2transcript[gene] = best_transcript

        if len(transcripts) > 1:
            print "##"
            print "## WARNING: Gene '%s' maps to %s transcripts: %s" % (gene, len(transcripts), transcripts)
            print "## WARNING: Choose transcript '%s' since it has the most mutation calls." % best_transcript

    return gene2transcript

def import_annotated_maf(infile):
    """Import the maf2maf (VEP) annotated MAF file.

    Notes
    -----
    (1) This infile is a tab-delimited MAF file. For full-spec, see:

            https://github.com/mskcc/vcf2maf/blob/master/docs/vep_maf_readme.txt

    (2) Replace 'Hugo_Symbol' to 'Gene_Symbol' (if not already done so).

    (3) This infile should additionally have the custom 'Source_Maf' column
        which indicate the (pre-merge) source MAF that each variant call
        originally came from.
    """

    annotated_maf_df = pd.read_table(infile, sep="\t", dtype=str, comment="#", header = 0)

    if 'Gene_Symbol' == list(annotated_maf_df)[0]:
        assert 'Hugo_Symbol' not in list(annotated_maf_df)
    else:
        assert 'Hugo_Symbol' == list(annotated_maf_df)[0]

        # Rename the 'Hugo_Symbol' column to 'Gene_Symbol' to avoid confusion. Reason
        # is that the symbol values in the 'Hugo_Symbol' can originate from various
        # gene_symbol sources and does not necessary need to come from 'HGNC'.
        # For details see:
        #    http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html
        #
        # The vep tool itself calls the field 'SYMBOL'. It is the vcf2maf tool
        # which maps these values to the 'Hugo_Symbol' field. See:
        #   https://github.com/mskcc/vcf2maf/blob/v1.6.9/vcf2maf.pl#L537
        annotated_maf_df = annotated_maf_df.rename(columns={'Hugo_Symbol': 'Gene_Symbol'})

    return annotated_maf_df

def import_transcriptcodon2trinuc_map(infile, keep_only_these_genes = None):
    """Import the transcript-codon to trinucleotides mapping into a 2-layer
    dictionary.

    Notes
    -----
    If user specified an input 'keep_only_these_genes' set, the only import the
    gene that appear in the set (this is done to reduce the memory footprint).
    """

    line_num = 0
    transcriptcodon2trinucleotides = dict()

    header_line = True

    expected_header_cols = ["Gene_Symbol", "Transcript_ID", "Codon_Pos", "Strand",
                            "Frame", "Chrom", "Position", "Trinucleotide"]

    for line in open(infile, 'rU'):

        cols = line.rstrip('\n').split('\t')

        if line.startswith("#"):
            continue

        if line.startswith("Gene_Symbol"):
            # Consistency check
            assert cols == expected_header_cols
            continue

        line_num += 1

        if line_num % 1000000 == 0:
            print "##"
            print "## PROGRESS: Processing transcriptcodon2trinuc_map line_num: %s." % line_num

        assert len(cols) == 8

        gene_symbol = cols[0]

        if keep_only_these_genes != None:
            if gene_symbol not in keep_only_these_genes:
                continue

        transcript_id = cols[1]
        codon_pos = int(cols[2])
        strand = cols[3]
        frame = int(cols[4])
        chrom = cols[5]
        position = int(cols[6])
        trinucleotide = cols[7]

        assert strand in ['-', '+']
        assert frame in [0,1,2]
        assert len(trinucleotide) == 3

        if transcript_id not in transcriptcodon2trinucleotides:
            transcriptcodon2trinucleotides[transcript_id] = dict()

        if codon_pos not in transcriptcodon2trinucleotides[transcript_id]:
            transcriptcodon2trinucleotides[transcript_id][codon_pos] = dict()

        assert frame not in transcriptcodon2trinucleotides[transcript_id][codon_pos]

        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame] = dict()
        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame]['Gene_Symbol'] = gene_symbol
        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame]['Strand'] = strand
        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame]['Chrom'] = chrom
        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame]['Position'] = position
        transcriptcodon2trinucleotides[transcript_id][codon_pos][frame]['Trinucleotide'] = trinucleotide

    # Consistency check
    for transcript in transcriptcodon2trinucleotides:

        protein_length = max(transcriptcodon2trinucleotides[transcript])

        assert set(range(1,protein_length+1)) == set(transcriptcodon2trinucleotides[transcript].keys())

    return transcriptcodon2trinucleotides

def complement_nt_letter(letter):
    """Get the complement of the input necleotide letter."""

    complement_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    return complement_map[letter]

def normalize_trinucleotide(trinucleotide):
    """Return the normalized representation of the input trinucleotide sequence

    Notes
    -----
    Each trinucleotide sequence has two possible representations (the sequence
    and its reverse complement). For example, 5'-ACG-3' and 5'-CGT-3' are two
    representations of the same trinucleotide sequence. To prevent ambiguity,
    choose the representation where the central nucleotide of the trinucleotide
    context is a C or a T is chosen.
    """

    # Consistency checks
    assert len(trinucleotide) == 3
    for letter in trinucleotide:
        assert letter in ['A', 'C', 'G', 'T']

    complement_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    reverse_complement = ""
    for letter in trinucleotide[::-1]:
        reverse_complement += complement_map[letter]

    # Consistency checks
    assert len(reverse_complement) == 3
    for letter in reverse_complement:
        assert letter in ['A', 'C', 'G', 'T']

    # Choose the seq where the middle nucleotide is a 'C' or a 'T'
    if trinucleotide[1] in ['C', 'T']:
        return trinucleotide
    elif reverse_complement[1] in ['C', 'T']:
        return reverse_complement
    else:
        raise Exception("Unexpected error.")

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
            return return_list

def get_hypermutated_samples(maf_df, source_maf_id2num_samples):
    """Identify the list of hypermutated_samples.

    Notes
    -----
    (1) For each panel (source_maf_id), construct the # mutations per sample distribution.
        Compute the Q1, median, Q3 and inter-quartile range (IQR) of the distribution.
        Define a hyptermutated sample as as 'right' outlier samples under the Tukey's
        method (# mutation > Q3 + 3x*QR).

    (2) The outlier detection is due seperately for each panel since each panel
        have a different target footprint (few hundred genes to exome to WGS). Hence
        distribution of mutations per sample can be vastly different between panels.

    (3) Also update source_maf_id2num_samples
    """

    hypermutated_samples = set()

    source_maf_ids = sorted(maf_df['Source_Maf'].unique())
    assert set(source_maf_id2num_samples.keys()) == set(source_maf_ids)

    for source_maf_id in source_maf_ids:

        num_sequenced_samples = source_maf_id2num_samples[source_maf_id]['Total']

        tmp_df = maf_df[maf_df['Source_Maf'] == source_maf_id]
        num_mutations_df = pd.DataFrame(tmp_df['Tumor_Sample_Barcode'].value_counts().rename('Num_Mutations'))

        num_maf_samples = len(num_mutations_df)

        # Add samples with no mutations to 'num_mutations_df'. This is necessary
        # to ensure that Q1, median, and Q3 are computed correctly.
        no_mutation_data_list = list()
        no_mutation_sample_list = list()
        for index in range(num_sequenced_samples - num_maf_samples):
            no_mutation_data = dict()
            no_mutation_data['Num_Mutations'] = 0
            no_mutation_sample = '%s_No_Mutation_Sample_%04d' % (source_maf_id, index + 1)
            no_mutation_data_list.append(no_mutation_data)
            no_mutation_sample_list.append(no_mutation_sample)
        no_mutations_df = pd.DataFrame(no_mutation_data_list, index=no_mutation_sample_list)
        num_mutations_df = pd.concat([num_mutations_df, no_mutations_df], axis=0)

        # Compute Statistics
        median = num_mutations_df['Num_Mutations'].median()
        q1_val = num_mutations_df['Num_Mutations'].quantile(q=0.25)
        q3_val = num_mutations_df['Num_Mutations'].quantile(q=0.75)
        interquartile_range = q3_val - q1_val
        outlier_cutoff = q3_val + 3 * interquartile_range
        hypermutated_df = num_mutations_df[num_mutations_df['Num_Mutations'] > outlier_cutoff]
        num_hypermutated_samples = len(hypermutated_df)

        assert hypermutated_df.index.is_unique

        verbose = False
        if verbose:
            mean = num_mutations_df['Num_Mutations'].mean()
            max_val = num_mutations_df['Num_Mutations'].max()

            frac_hypermutated_samples = num_hypermutated_samples / float(num_sequenced_samples)
            num_total_mutations = num_mutations_df['Num_Mutations'].sum()
            num_hypermutated_mutations = hypermutated_df['Num_Mutations'].sum()
            frac_hypermutated_mutations = num_hypermutated_mutations / float(num_total_mutations)
            print "## -----------------------------------------------------------------"
            print "## Get HyperMutated Samples Summary (source_maf_id: '%s')" % source_maf_id
            print "##   q1_val: %5.1f" % q1_val
            print "##   median: %5.1f" % median
            print "##   mean: %5.1f" % mean
            print "##   q3_val: %5.1f" % q3_val
            print "##   max_val: %5.1f" % max_val
            print "##   outlier_cutoff : %5.1f" % outlier_cutoff
            print "##   num_total_samples (MAF): %s" % num_maf_samples
            print "##   num_total_samples (Sequenced): %s" % num_sequenced_samples
            print "##   num_hypermutated_samples: %s/%s (%.2f%%)" % (num_hypermutated_samples,
                                                                     num_sequenced_samples,
                                                                     100*frac_hypermutated_samples)
            print "##   num_hypermutated_mutations: %s/%s (%.2f%%)" % (num_hypermutated_mutations,
                                                                       num_total_mutations,
                                                                       100*frac_hypermutated_mutations)
            print "## -----------------------------------------------------------------"

        # Update source_maf_id2num_samples
        source_maf_id2num_samples[source_maf_id]['Hypermutated'] = num_hypermutated_samples
        source_maf_id2num_samples[source_maf_id]['NotHypermutated'] = num_sequenced_samples - num_hypermutated_samples

        perpanel_hypermutated_samples = set(hypermutated_df.index.values)

        # Consistency that there is no duplicated sample_id across panels.
        assert len(hypermutated_samples & perpanel_hypermutated_samples) == 0
        hypermutated_samples.update(perpanel_hypermutated_samples)

    return hypermutated_samples


def compute_q_values_wrapper(out_data_list, total_enumerated_codons):
    """
    Computing qvalues for the Logistic Regression Hotspot Detection tests and
    add the compute values to the inputted out_data_list

    Notes
    -----
    (1) Pad the pvalues vector with 1.0 for codons that were enumerated BUT not
        evaluated (fewer mutation calls than 'min_count',
    """

    print "##"
    print "## Computing Q-Values for Logistic Regression Hotspot Detection Test..."

    assert total_enumerated_codons >= len(out_data_list)

    pvalues = list()
    for idx in range(total_enumerated_codons):

        if idx < len(out_data_list):
            pvalue = out_data_list[idx]['PValue']
        else:
            pvalue = 1.0

        pvalues.append(pvalue)

    qvalues = compute_q_values(pvalues)

    # Consistency checkes
    assert total_enumerated_codons == len(pvalues)
    assert total_enumerated_codons == len(qvalues)

    for idx in range(len(out_data_list)):
        out_data_list[idx]['QValue'] = qvalues[idx]

    return None


def compute_q_values(pvalues):
    """
    Compute the FDR-adjusted p-values (q-values) from the input list of
    uncorrected p-values (floats).

    Return the q-values as a list of floats.

    Note
    ----
    (1) Use the Benjamini  and Hochberg (1995) method, which assumes independence
        or positive correlation between the p-values.

    (2) Take care to maintain the order of the elements in the list.

    (3) Tested that this function can correctly handle the edge case where one or
        more pvalue(s) equal 0.0.
    """

    # Consistency check.
    for idx in range(len(pvalues)):

        pvalue = pvalues[idx]

        assert isinstance(pvalue, float)
        assert pvalue >= 0.0
        assert pvalue <= 1.0

    pvalues = np.asarray(pvalues)

    # fdr_bh: Benjamini/Hochberg  (non-negative)
    _, qvalues, _, _ = statsmodels.sandbox.stats.multicomp.multipletests(pvalues,
                                                                         method='fdr_bh',
                                                                         is_sorted=False)

    qvalues = list(qvalues)

    # Consistency check
    for qvalue in qvalues:
        assert isinstance(qvalue, float)
        assert qvalue >= 0.0
        assert qvalue <= 1.0

    return qvalues

def map_source_maf_id_to_hardcoded_num_samples(source_maf_ids, maf_df):
    """Create a dictionary mapping source_maf_id to # of unique samples with the
    source_maf_id in the maf_df.

    Note
    ----
    (1) This function is meant only for local testing purposes. Don't use to
        actual production analysis.

    (2) Hard-coded numbers which also counts samples with no mutations (and hence
        do not how up in the MAF file).

    (3) Updated on April 27th, 2017. Inferred from # samples in the Clinical Data
        File (private/merged.clinical_data.r13.tsv) [EXCLUDE Male + Metastatic]
    """

    no_mutation_DFCI_all = 18  # Num samples that were sequenced but no mutation called (all 3 panels)
    no_mutation_DFCI_pieces = divide_num_by_weights_and_round(no_mutation_DFCI_all, [127, 548, 75])

    no_mutation_MSK_all = 15  # Num samples that were sequenced but no mutation called (all 2 panels)
    no_mutation_MSK_pieces = divide_num_by_weights_and_round(no_mutation_MSK_all, [155, 161])

    no_mutation_METABRIC = 138  # Num samples that were sequenced but no mutation called (1 panel)

    source_maf2value = dict()
    source_maf2value['dfci.DFCI-ONCOPANEL-1.r6.vep.maf'] = 203 + no_mutation_DFCI_pieces[0]
    source_maf2value['dfci.DFCI-ONCOPANEL-2.r6.vep.maf'] = 858 + no_mutation_DFCI_pieces[1]
    source_maf2value['dfci.DFCI-ONCOPANEL-3.r6.vep.maf'] = 248 + no_mutation_DFCI_pieces[2]
    # source_maf2value['genie_msk.IMPACT341.r3.vep.maf'] = 155 +  no_mutation_MSK_pieces[0]
    # source_maf2value['genie_msk.IMPACT410.r3.vep.maf'] = 161 +  no_mutation_MSK_pieces[1]
    # source_maf2value['metabric.vep.maf'] = 2368 + no_mutation_METABRIC
    # source_maf2value['sanger_somatic_coding_annotated.maf'] = 548 + 0
    # source_maf2value['tcga.cell_2015.r2.vep.maf'] = 808 + 0

    keys = source_maf2value.keys()

    # Get num_samples, including both sample WITH AND WITHOUT mutations.
    source_maf_id2num_samples = dict()
    for source_maf_id in source_maf2value:
        if source_maf_id in source_maf_ids:
            source_maf_id2num_samples[source_maf_id] = dict()
            source_maf_id2num_samples[source_maf_id]['Total'] = source_maf2value[source_maf_id]

    # Consistency checks
    assert set(source_maf_id2num_samples.keys()) == set(source_maf_ids)
    source_maf_id2num_maf_samples = map_source_maf_id_to_num_maf_samples(source_maf_ids, maf_df)
    for source_maf_id in source_maf_id2num_samples:
        assert source_maf_id2num_samples[source_maf_id]['Total'] >= source_maf_id2num_maf_samples[source_maf_id]['Total']

    return source_maf_id2num_samples

def divide_num_by_weights_and_round(in_num, weights):
    """Divide up the input value (an integer) by the pieces according to the
    specified weights and round the return values.

    Notes
    -----
    (1) Round the return values to whole number in such a way that:
         (a) Ensure that the rounded_vals still sum up to in_num
         (b) Minimize the absolute deviation between rounded_vals and unrounded_vals

    (2) The implemention below search for the solution that minimize the
        absolute deviation via brute force. It gets the job done but is certainly
        not the most elegant + optimized solution.
    """

    assert isinstance(in_num, int)
    for weight in weights:
        assert isinstance(weight, int) or isinstance(weight, float)

    num_vals = len(weights)

    # Normalized the weights
    norm_weights = [weight/float(sum(weights)) for weight in weights]

    unrounded_vals = [in_num * norm_weight for norm_weight in norm_weights]


    initial_guess_rounded_vals = [int(val) for val in unrounded_vals]

    max_search_diff = 2
    search_diff_range = range(-max_search_diff, max_search_diff + 1)

    best_rounded_vals = None
    best_total_abs_deviation = None

    for diffs in itertools.product(search_diff_range, repeat=num_vals):
        rounded_vals = initial_guess_rounded_vals[:]  # Hard-copy
        total_abs_deviation = 0
        for index in range(num_vals):
            rounded_vals[index] += diffs[index]
            total_abs_deviation += abs(rounded_vals[index] - unrounded_vals[index])

        if sum(rounded_vals) == in_num:
            if best_rounded_vals == None or total_abs_deviation < best_total_abs_deviation:
                best_rounded_vals = rounded_vals
                best_total_abs_deviation = total_abs_deviation

    # Consistent checks
    assert best_rounded_vals != None
    assert best_total_abs_deviation != None
    assert num_vals == len(best_rounded_vals)
    for val in best_rounded_vals:
        assert isinstance(val, int)

    return best_rounded_vals

#######################################################################################
################ Functions That Originated From Mutation Counts Script ################
#######################################################################################

def map_source_maf_id_to_num_maf_samples(source_maf_ids, maf_df):
    """Create a dictionary mapping source_maf_id to # of unique samples with the
    source_maf_id in the maf_df."""

    return_dict = dict()

    # Consistency check
    assert isinstance(source_maf_ids, list)
    assert len(source_maf_ids) == len(set(source_maf_ids))

    for source_maf_id in source_maf_ids:
        samples_in_maf = get_unique_samples(source_maf_id, maf_df)
        return_dict[source_maf_id] = dict()
        return_dict[source_maf_id]['Total'] = len(samples_in_maf)

    # Print a WARNING message if there are orphan 'Source_Maf' values in
    # maf_df that is not in source_maf_ids.
    orphan_maf_ids = list(set(maf_df['Source_Maf'].unique()) - set(source_maf_ids))
    if len(orphan_maf_ids) != 0:
        print "##"
        print "## WARNING: Found orphan 'Source_Maf' value(s) %s which do not appear" % orphan_maf_ids,
        print "in provided source_maf_ids list."
        print "## WARNING: Samples in MAF file with these 'Source_Maf' value(s)",
        print "will be excluded from the targeted/untargeted sample counts."

    return return_dict

def get_unique_samples(source_maf_id, maf_df):
    """Get a set of unique samples in maf_df with the specified source_maf_id."""

    temp_df = maf_df[maf_df['Source_Maf'] == source_maf_id]

    samples = set(temp_df['Tumor_Sample_Barcode'].unique())

    if len(samples) == 0:
        print "##"
        print "## WARNING: Found 0 sample in MAF file with source_maf_id %s." % repr(source_maf_id)

    return samples

def map_source_maf_id_to_center(source_maf_ids, maf_df):
    """Create a dictionary mapping source_maf_id to center name.

    Note
    ----
    Check that if two rows in maf_df has the same 'source_maf_id' value
    then the two rows must have the same 'center' value as well.
    """

    return_dict = dict()

    # Consistency check
    assert isinstance(source_maf_ids, list)
    assert len(source_maf_ids) == len(set(source_maf_ids))

    for source_maf_id in source_maf_ids:
        temp_df = maf_df[maf_df['Source_Maf'] == source_maf_id]
        centers = list(temp_df['Center'].unique())

        if len(centers) == 0:
            raise Exception("source_maf_id '%s' maps to no center." % source_maf_id)

        if len(centers) > 1:
            raise Exception("source_maf_id '%s' maps to multiple centers: %s" % (source_maf_id, centers))

        return_dict[source_maf_id] = centers[0]

    return return_dict

def apply_missense_snp_filter(maf_df):
    """Filter for SNP variants that have a 'missense' effect.

    Notes
    -----
    (1) Previously applied the following filter:

            maf_df = maf_df[maf_df['Variant_Classification'] == 'Missense_Mutation']

        However found that the above filter was too permissive and include (rare)
        variants where 'Consequence' value is 'coding_sequence_variant' (which is
        not a 'missense_variant').

        The reason is because the 'Missense_Mutation' value is computed by vcf2maf
        tool and maps one of the following 'Consequence' values:
            (1) missense_variant
            (2) coding_sequence_variant
            (3) conservative_missense_variant
            (4) rare_amino_acid_variant

         Also See:
            http://useast.ensembl.org/info/genome/variation/predicted_data.html
            http://www.sequenceontology.org/miso/current_release/term/SO:0001580
    """

    maf_df = maf_df[maf_df['Consequence'].str.contains('missense_variant')]
    maf_df = maf_df[maf_df['Variant_Type'] == 'SNP']

    return maf_df

def apply_exac_filter(maf_df):
    """Filter-out variant calls with Adjusted ExAC Allele Frequency above 0.06% (to remove
    putative germline variants)."""

    total_variant_calls = len(maf_df)

    # Print excluded variants for debugging/logging purpose.
    exclude_indices = maf_df[maf_df['ExAC_AF_Adj'].astype(float).fillna(0.0) > 0.0006].index.tolist()

    print "##"
    print "## Performing ExAC Allele Frequency Filter..."
    print "##"
    print "## Filtered-out %s variant calls with ExAC_AF_Adj > 0.06%%" % format(len(exclude_indices), ',d'),
    print "(Note: out of %s total variant calls)" % format(total_variant_calls, ',d')

    # Drop the variants
    maf_df.drop(exclude_indices, inplace=True)

    # Consistency check
    assert len(maf_df) == total_variant_calls - len(exclude_indices)

    return maf_df

######################################################
################ Target BED Functions ################
######################################################

def map_source_maf_id_to_target_intervals(source_maf_ids, target_bed_files):
    """Create a dictionary mapping source_maf_id to target_bed_file. Assumes that
    the two lists are ordered such that their index can be used for the mapping."""

    return_dict = dict()

    # Consistency checks
    assert isinstance(source_maf_ids, list)
    assert isinstance(target_bed_files, list)
    assert len(source_maf_ids) == len(set(source_maf_ids))
    assert len(source_maf_ids) == len(target_bed_files)

    for index in range(len(source_maf_ids)):
        source_maf_id = source_maf_ids[index]
        target_intervals = import_target_intervals(target_bed_files[index])
        return_dict[source_maf_id] = target_intervals

    return return_dict

def import_target_intervals(target_bed_file):
    """Import the target_intervals from the target_bed_file."""

    bed_obj = pybedtools.BedTool(target_bed_file)
    target_intervals = bed_obj.as_intervalfile()

    return target_intervals

def create_point_interval(chrom, coordinate):
    """Create pybedtools interval for a point genomics position."""

    interval = pybedtools.create_interval_from_list([chrom, coordinate, coordinate])

    return interval

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to input (multi-center) MAF file.")

    parser.add_argument("--in_map", action="store", required=True,
                        metavar='FILE',
                        help="Path to input codon2trinucleotide map file")

    parser.add_argument("--in_source_maf_ids", action="store", required=True,
                        nargs='+',
                        help="Source_Maf IDs for mapping with target bed files (by ordering).")

    parser.add_argument("--in_target_beds", action="store", required=True,
                        nargs='+',
                        help="Target bed files for mapping with Source_Maf IDs (by ordering).")

    parser.add_argument("--min_count", action="store", type=int, default=3,
                        help="Don't evaluate codon with mutation counts below this threshold.")

    parser.add_argument("--use_toy_stats_model", action="store_true", default=False,
                        help="Use Simplified Logistic Regression Model (i.e. No Fixed Effects).")

    parser.add_argument("--out_table", action="store", required=True,
                        metavar='FILE',
                        help="""Path to the output table.""")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_maf:", repr(options.in_maf)
    print "##   in_map:", repr(options.in_map)
    print "##   in_source_maf_ids:"
    for index, source_maf_id in enumerate(options.in_source_maf_ids, start=1):
        print "##       (%d) %s" % (index, repr(source_maf_id))
    print "##   in_target_beds:"
    for index, target_bed in enumerate(options.in_target_beds, start=1):
        print "##       (%d) %s" % (index, repr(target_bed))
    print "##   min_count:", repr(options.min_count)
    print "##   use_toy_stats_model:", repr(options.use_toy_stats_model)
    print "##   out_table:", repr(options.out_table)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
