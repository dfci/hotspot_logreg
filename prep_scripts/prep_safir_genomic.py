#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the SAFIR data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

   - Study is a whole-exome mutational profiling of metastatic breast cancer.
        - hg19 build
        - strand not specified. I'm assuming Positive

   - Paper: https://www.ncbi.nlm.nih.gov/pubmed/28027327
   - Public dataset downloaded from cBioPortal datahub.

EXAMPLES

    ./prep_safir_genomic.py \
        --in-maf ./data_mutations_extended.txt \
        --out-genomic ./safir.maf

AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (May-2017)
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


class Genomic:

    def __init__(self, args):

        self.args = args
        self.genomic_df = None
        self.maf_df = pd.DataFrame(columns=MAF_columns)

        self.num_rows = None
        self.num_samples = None
        self.i = 0

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.genomic_df = pd.read_csv(self.args.in_genomic, sep='\t')
        self.num_rows = len(self.genomic_df.index)
        self.num_samples = len(self.genomic_df.Tumor_Sample_Barcode.unique())

    def write_results(self):
        """Write clinical dataframe to file."""

        # clinical file
        with open(self.args.out_genomic, 'w') as ff:
            ff.write('#version 1.0\n')
            ff.write('#Original .maf contained %d rows and %d samples\n' % (self.num_rows, self.num_samples))
            self.maf_df.to_csv(ff, sep='\t', index=False)

    def create_maf(self):

        print '## Input contains %d rows\n' % self.num_rows
        self.maf_df.Tumor_Sample_Barcode = self.genomic_df.Tumor_Sample_Barcode
        self.maf_df.Center = 'SAFIR'
        self.maf_df.NCBI_Build = 'GRCh37'
        self.maf_df.Chromosome = self.genomic_df.Chromosome
        self.maf_df.Start_Position = self.genomic_df.Start_Position
        self.maf_df.End_Position = self.genomic_df.End_Position
        self.maf_df.Strand = '+'
        self.maf_df.Reference_Allele = self.genomic_df.Reference_Allele
        self.maf_df.Tumor_Seq_Allele1 = self.genomic_df.Tumor_Seq_Allele1
        self.maf_df.Tumor_Seq_Allele2 = self.genomic_df.Tumor_Seq_Allele2
        self.maf_df = self.maf_df.fillna(value='', axis=1)


def main(args):

    run = Genomic(args)
    run.load_data()
    run.create_maf()
    run.write_results()



if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-maf', dest='in_genomic', required=True,
                        help='Path to input genomic .maf file from cBioPortal')

    parser.add_argument('--out-genomic', dest='out_genomic', required=True,
                        help='Path to output genomic file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
