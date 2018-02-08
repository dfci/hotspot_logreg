#!/usr/bin/env python2

"""
SYNOPSIS

    For each Panel (Source_Maf), create a '# mutations per sample' histogram and
    determine an outlier cutoff value to detect hypermutated samples.

NOTES

    (1) Output the plot for each panel to 'hypermutated_samples.<source_maf_id>.png'

    (2) Compute the Q1, median, Q3 and inter-quartile range (IQR) of the distribution.
        Define a hyptermutated sample as as 'right' outlier samples under the Tukey's
        method (# mutation > Q3 + 3 x IQR).

    (3) The outlier detection is due seperately for each panel since each panel
        have a different target footprint (few hundred genes to exome to WGS). Hence
        distribution of mutations per sample can be vastly different between panels.

EXAMPLES

    ./hypermutated_samples.py \
        --in_maf ../../private/merged.dfci_plus_public.r6.filtered_clinical.vep.maf

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import argparse
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
import math
import itertools
import matplotlib.pyplot as plt

def main(options):

    maf_df = import_annotated_maf(options.in_maf)

    maf_df = apply_missense_snp_filter(maf_df)

    maf_df = apply_exac_filter(maf_df)

    source_maf_id2num_samples = get_source_maf_id_to_num_samples()

    for source_maf in maf_df['Source_Maf'].unique():
        run(options, maf_df, source_maf, source_maf_id2num_samples)

def run(options, maf_df, source_maf, source_maf_id2num_samples):

    plot_label = get_label(source_maf)

    num_sequenced_samples = source_maf_id2num_samples[source_maf]

    maf_df = maf_df[maf_df['Source_Maf'] == source_maf]

    num_muts_df = pd.DataFrame(maf_df['Tumor_Sample_Barcode'].value_counts().rename('Num_Mutations'))

    # Quick Hack/Check
    if plot_label == 'TCGA-BRCA':
        num_muts_df = num_muts_df[num_muts_df['Num_Mutations'] < 3000]

    num_maf_samples = len(num_muts_df)

    assert num_sequenced_samples >= num_maf_samples

    # Add Sequenced Samples with no mutations to 'num_muts_df'
    no_mutation_data_list = list()
    no_mutation_sample_list = list()
    for index in range(num_sequenced_samples - num_maf_samples):
        no_mutation_data = dict()
        no_mutation_data['Num_Mutations'] = 0
        no_mutation_sample = 'No_Mutation_Sample_%04d' % (index + 1)
        no_mutation_data_list.append(no_mutation_data)
        no_mutation_sample_list.append(no_mutation_sample)
    no_muts_df = pd.DataFrame(no_mutation_data_list, index=no_mutation_sample_list)
    num_muts_df = pd.concat([num_muts_df, no_muts_df], axis=0)

    mean = num_muts_df['Num_Mutations'].mean()
    median = num_muts_df['Num_Mutations'].median()
    q1_val = num_muts_df['Num_Mutations'].quantile(q=0.25)
    q3_val = num_muts_df['Num_Mutations'].quantile(q=0.75)
    interquartile_range = q3_val - q1_val
    max_val =  num_muts_df['Num_Mutations'].max()
    outlier_cutoff = q3_val + 3 * interquartile_range
    outlier_num_muts_df = num_muts_df[num_muts_df['Num_Mutations'] > outlier_cutoff]

    num_hypermutated_samples = len(outlier_num_muts_df)
    frac_hypermutated_samples = num_hypermutated_samples / float(num_sequenced_samples)

    num_total_mutations = num_muts_df['Num_Mutations'].sum()
    num_hypermutated_mutations = outlier_num_muts_df['Num_Mutations'].sum()
    frac_hypermutated_mutations = num_hypermutated_mutations / float(num_total_mutations)

    num_bins = min(70, int(max_val))

    count_bins, bin_edges = plot_histrogram(plt, num_muts_df, num_bins,
                                            q1_val, median, q3_val, outlier_cutoff)

    xmin = min(0.0, -0.5 * bin_edges[1])
    xmax = max(max_val * 1.1, outlier_cutoff * 5.00)
    # xmax = min(max_val * 1.1, outlier_cutoff * 3.00)
    ymin = 0
    ycutoff = 40
    ymax = max(max(count_bins) * 1.1, ycutoff)

    print "## -----------------------------------------------------------------"
    print "## Summary Stats (source_maf %s)" % source_maf
    print "##   num_bins: %d | xmin: %d | xmax: %d | ymin: %d | ycutoff: %d | ymax: %d" % (num_bins, xmin, xmax, ymin, ycutoff, ymax)
    print "##   q1_val: %5.1f" % q1_val
    print "##   median: %5.1f" % median
    print "##   mean: %5.1f" % mean
    print "##   q3_val: %5.1f" % q3_val
    print "##   max_val: %5.1f" % max_val
    print "##   outlier_cutoff : %5.1f" % outlier_cutoff
    print "##   num_total_samples (MAF): %s" % num_maf_samples
    print "##   num_total_samples (Hard-coded): %s" % num_sequenced_samples
    print "##   num_hypermutated_samples: %s/%s (%.2f%%)" % (num_hypermutated_samples, num_sequenced_samples, 100*frac_hypermutated_samples)
    print "##   num_hypermutated_mutations: %s/%s (%.2f%%)" % (num_hypermutated_mutations, num_total_mutations, 100*frac_hypermutated_mutations)
    print "## -----------------------------------------------------------------"


    plt.clf()


    plt.rcParams['figure.figsize'] = (10, 10)
    params = {'legend.fontsize': 12, 'axes.labelsize':16, 'axes.labelcolor':  "black",
              'axes.labelweight': "bold", "xtick.labelsize": 16, "ytick.labelsize": 16,
              'legend.markerscale': 1.5}
    plt.rcParams.update(params)

    if max(count_bins) < 50:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1
        count_bins, bin_edges = plot_histrogram(ax1, num_muts_df, num_bins,
                                                q1_val, median, q3_val, outlier_cutoff)

        plt.ylim(ymin, ymax)

    else:

        # If we were to simply plot pts, we'd lose most of the interesting
        # details due to the outliers. So let's 'break' or 'cut-out' the y-axis
        # into two portions - use the top (ax) for the outliers, and the bottom
        # (ax2) for the details of the majority of our data
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
        gs.update(hspace=0.05)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        # plot the same data on both axes
        for curr_ax in [ax1, ax2]:
            count_bins, bin_edges = plot_histrogram(curr_ax, num_muts_df, num_bins,
                                                    q1_val, median, q3_val, outlier_cutoff)

        ax1.set_ylim(ycutoff, ymax)  # High-counts
        ax2.set_ylim(ymin, ycutoff)  # Low-counts

        # hide the spines between ax and ax2
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.xaxis.tick_top()
        ax1.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()

        # Add Diagonal
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    ax1.set_xlim(xmin, xmax)
    ax2.set_xlim(xmin, xmax)
    ax1.legend(loc="upper right")
    ax1.set_title('Number of Missense Mutations per Sample', fontsize = 20, fontweight='bold')
    ax2.set_xlabel('Number of Missense Mutations')
    fig.text(0.06, 0.5, 'Number of Samples', ha='center', va='center', rotation='vertical', fontsize = 16, weight= 'bold')

    # plt.tight_layout()
    plt.savefig('hypermutated_samples.%s.png' % plot_label)
    plt.clf()

def plot_histrogram(matplotlib_obj, num_muts_df, num_bins, q1_val, median, q3_val, outlier_cutoff):

    count_bins, bin_edges, _ = matplotlib_obj.hist(num_muts_df['Num_Mutations'], num_bins, facecolor='green', alpha=0.75, label=None, align='left')


    matplotlib_obj.axvline(x=q1_val, ymin=0, ymax = max(count_bins), linewidth=2, color='b', linestyle='dashed', label='Q1: %5.1f' % q1_val)
    matplotlib_obj.axvline(x=median, ymin=0, ymax = max(count_bins), linewidth=2, color='b', label='Median: %5.1f' % median)
    matplotlib_obj.axvline(x=q3_val, ymin=0, ymax = max(count_bins), linewidth=2, color='b', linestyle='dashed', label='Q3: %5.1f' % q3_val)
    matplotlib_obj.axvline(x=outlier_cutoff, ymin=0, ymax = max(count_bins), linewidth=2, color='r', label='Q3 + 3*IQR: %5.1f' % outlier_cutoff)

    # outlier_num_muts_df = num_muts_df[num_muts_df['Num_Mutations'] > outlier_cutoff]
    # matplotlib_obj.hist(outlier_num_muts_df['Num_Mutations'], bin_edges, facecolor='red', alpha=0.75)

    return count_bins, bin_edges

def get_label(source_maf):

    source_maf2label = dict()
    source_maf2label['dfci.DFCI-ONCOPANEL-1.r4.vep.maf'] = "DFCI_ONCOPANEL_v1"
    source_maf2label['dfci.DFCI-ONCOPANEL-2.r4.vep.maf'] = "DFCI_ONCOPANEL_v2"
    source_maf2label['dfci.DFCI-ONCOPANEL-3.r4.vep.maf'] = "DFCI_ONCOPANEL_v3"
    source_maf2label['genie_msk.IMPACT341.r3.vep.maf'] = "MSK_IMPACT341"
    source_maf2label['genie_msk.IMPACT410.r3.vep.maf'] = "MSK_IMPACT410"
    source_maf2label['metabric.vep.maf'] = "METABRIC"
    source_maf2label['sanger_somatic_coding_annotated.maf'] = "SangerWGS"
    source_maf2label['tcga.cell_2015.r2.vep.maf'] = "TCGA-BRCA"
    label = source_maf2label[source_maf]

    return label

def get_source_maf_id_to_num_samples():
    """Create and return dictionary mapping # of Sequenced Samples for each Source_Maf (Panel).

    Notes
    -----
      (1) This includes the # of Samples with no mutations (and hence does not
          show up in the MAF file.

      (2) Updated on April 27th, 2017. Inferred from # samples in the Clinical Data File.
    """

    no_mutation_DFCI_all = 18  # Num samples that were sequenced but no mutation called (all 3 panels)
    no_mutation_DFCI_pieces = divide_num_by_weights_and_round(no_mutation_DFCI_all, [127, 548, 75])

    no_mutation_MSK_all = 15  # Num samples that were sequenced but no mutation called (all 2 panels)
    no_mutation_MSK_pieces = divide_num_by_weights_and_round(no_mutation_MSK_all, [155, 161])

    no_mutation_METABRIC = 138  # Num samples that were sequenced but no mutation called (1 panel)

    source_maf2value = dict()
    source_maf2value['dfci.DFCI-ONCOPANEL-1.r4.vep.maf'] = 127 + no_mutation_DFCI_pieces[0]
    source_maf2value['dfci.DFCI-ONCOPANEL-2.r4.vep.maf'] = 548 + no_mutation_DFCI_pieces[1]
    source_maf2value['dfci.DFCI-ONCOPANEL-3.r4.vep.maf'] = 75 + no_mutation_DFCI_pieces[2]
    source_maf2value['genie_msk.IMPACT341.r3.vep.maf'] = 155 +  no_mutation_MSK_pieces[0]
    source_maf2value['genie_msk.IMPACT410.r3.vep.maf'] = 161 +  no_mutation_MSK_pieces[1]
    source_maf2value['metabric.vep.maf'] = 2368 + no_mutation_METABRIC
    source_maf2value['sanger_somatic_coding_annotated.maf'] = 548 + 0
    source_maf2value['tcga.cell_2015.r2.vep.maf'] = 808 + 0

    print "##"
    print "## DEBUG: no_mutation_DFCI_pieces: %s" % no_mutation_DFCI_pieces
    print "## DEBUG: no_mutation_MSK_pieces: %s" % no_mutation_MSK_pieces
    print "## DEBUG: source_maf_id_to_num_samples:"
    for source_maf_id in sorted(source_maf2value.keys()):
        print "##      %s --> %s" % (source_maf_id, source_maf2value[source_maf_id])

    return source_maf2value

def apply_exac_filter(maf_df):
    """Filter-out variant calls with Adjusted ExAC Allele Frequency above 0.06%
    (to remove putative germline variants)."""

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

def import_annotated_maf(infile):
    """Import the maf2maf (VEP) annotated MAF file.

    Notes
    -----
    (1) This infile is a tab-delimited MAF file. For full-spec, see:

            https://github.com/mskcc/vcf2maf/blob/master/docs/vep_maf_readme.txt

    (2) This infile should additionally have the custom 'Source_Maf' column
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

if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to input (multi-center) MAF file.")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_maf:", repr(options.in_maf)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
