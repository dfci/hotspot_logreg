#!/usr/bin/env python2

"""
SYNOPSIS
    As of June 12th, 2017, Foundation Medicine has not released their
    proprietary .bed files. The Foundation medicine team has said that
    they fully tile all exons in a given gene list. This list is given in
    the file T5a_genes.txt.

    To create an approximation of the Foundation .bed file I began with the .bed
    file from the GENIE dataset downloaded here:
        https://www.synapse.org/#!Synapse:syn7844527

    I then took all exonic regions with a SEQ_ASSAY_ID of 'VICC-01-T5a' for all genes
    in the given gene list 'T5a_genes.txt'. There were a few genes that did not match between
    GENIE's .bed file and the given Foundation genelist because different synonyms of the same
    gene were used. These are converted in the GENIE dataframe to the Foundation values through
    the ParseBed.synonym_dict object.

EXAMPLES

    ./prep_foundation_bed.py \
        --in-tcga-bed ./genie_combined.bed \
        --in-genelist ./T5a_genes.txt \
        --out-bed ./foundation.v1.bed

AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (June-2017)
"""

import sys
import os
import json
import time
import argparse
import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)


class ParseBed:
    def __init__(self, args):
        self.args = args
        self.genelist = None
        self.genie_bed_df = None
        self.foundation_bed_df = None
        self.synonym_dict = {
            'PAK5': 'PAK7',
            'ADGRA2': 'GPR124',
            'MRE11': 'MRE11A',
            'EMSY': 'C11orf30'
        }
        self.skipped_genes = []

    def load_bed(self):
        self.genie_bed_df = pd.read_csv(self.args.in_bed, sep='\t')
        self.foundation_bed_df = pd.DataFrame(columns=self.genie_bed_df.columns)
        for gene in self.synonym_dict:
            self.genie_bed_df.replace(to_replace=gene, value=self.synonym_dict[gene], inplace=True)

    def load_genelist(self):
        df = pd.read_csv(self.args.in_genelist, sep='\t', header=-1)
        df.columns = ['gene_names']
        self.genelist = df.gene_names.unique().tolist()

    def subset_by_genelist(self):

        f_panel = (self.genie_bed_df.SEQ_ASSAY_ID == 'VICC-01-T5a')
        f_exon = (self.genie_bed_df.Feature_Type == 'exon')

        for gene in self.genelist:
            f1 = (self.genie_bed_df.Hugo_Symbol == gene)
            df = self.genie_bed_df[f1 & f_panel & f_exon]

            if df.size == 0:
                self.skipped_genes.append(gene)
                print '##\n## WARNING: %s was not found in the GENIE dataframe. It was excluded from' \
                      'the final .bed file.'
                continue

            self.foundation_bed_df = self.foundation_bed_df.append(df)

    def write_results(self):

        with open(self.args.out_bed, 'w') as ff:
            ff.write('#version 1.0\n#')
            self.foundation_bed_df.to_csv(ff, sep='\t', index=False)

        if len(self.skipped_genes) > 0:
            print '##\n## WARNING: The following genes were excluded from output .bed file:\n' \
                  '## %s\n' % ''.join(self.skipped_genes)


def main(args):
    run = ParseBed(args)
    run.load_bed()
    run.load_genelist()
    run.subset_by_genelist()
    run.write_results()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-bed', dest='in_bed', required=True,
                        help='Path to GENIE .bed file')

    parser.add_argument('--in-genelist', dest='in_genelist', required=True,
                        help='Path to Foundation gene list file')

    parser.add_argument('--out-bed', dest='out_bed', required=True,
                        help='Path to output Foundation .bed file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
