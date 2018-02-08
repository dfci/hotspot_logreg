#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the GENIE data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

   - Public dataset downloaded from Synapse:
    https://www.synapse.org/#!Synapse:syn7851246

EXAMPLES

    ./prep_genie_genomic.py \
        --in-maf ./data_mutations_extended.txt \
        --in-clinical ./ \
        --out-genomic ./

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

maf_columns = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
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

CENTERS = ['grcc', 'vicc', 'mda', 'uhn']


class Genomic:

    def __init__(self, args):

        self.args = args

        self.init_genomic_df = None

        template_dict = {
            'clinical_df': None,
            'init_genomic_df': None,
            'num_rows': None,
            'num_samples': None
        }
        self.genie_dict = {}
        for center in CENTERS:
            self.genie_dict[center] = template_dict.copy()

        self.i = 0

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.init_genomic_df = pd.read_csv(self.args.in_genomic, sep='\t', skiprows=1, low_memory=False)
        for center in CENTERS:

            filename = '%s/%s.clinical.r1.tsv' % (self.args.in_clinical, center.upper())
            clinical_df = pd.read_csv(filename, sep='\t', skiprows=1)

            self.genie_dict[center]['clinical_df'] = clinical_df
            self.genie_dict[center]['num_rows'] = len(clinical_df.index)
            self.genie_dict[center]['num_samples'] = len(clinical_df.Tumor_Sample_Barcode.unique())

    def write_results(self, df, center):
        """Write clinical dataframe to file."""

        num_rows = self.genie_dict[center]['num_rows']
        num_samples = self.genie_dict[center]['num_rows']
        filename = '%s/GENIE.%s.r1.maf' % (self.args.out_genomic, center)

        with open(filename, 'w') as ff:
            ff.write('#version 1.0\n')
            ff.write('#Original .maf contained %d rows and %d samples\n' % (num_rows, num_samples))
            df.to_csv(ff, sep='\t', index=False)

    def remove_non_breast_cancer_samples(self):
        """
        Using the output of prep_msk_clinical.py as a reference, remove all samples that did not have
        breast cancer oncotree codes.
        """

        for center in CENTERS:
            clinical_df = self.genie_dict[center]['clinical_df']

            keep_these_sample_ids = clinical_df.Tumor_Sample_Barcode.unique().tolist()

            f1 = (self.init_genomic_df.Tumor_Sample_Barcode.isin(keep_these_sample_ids))
            self.genie_dict[center]['init_genomic_df'] = self.init_genomic_df[f1]

    def remove_indels(self, df):
        """
        As of 05-19-2017, the decision was made to drop all indels from downstream analyses.
        This function does this by removing all rows where the lengths of the Reference_Allele,
        Tumor_Seq_Allele1, and Tumor_Seq_Allele2 are not equal. Note that this allows for any number
        of nucleotide polymorphisms, not just SNPs.
        """

        f1 = (len(df.Reference_Allele) == len(df.Tumor_Seq_Allele1))
        f2 = (len(df.Reference_Allele) == len(df.Tumor_Seq_Allele2))
        f3 = (df.Reference_Allele != '-')
        f4 = (df.Tumor_Seq_Allele1 != '-')
        f5 = (df.Tumor_Seq_Allele2 != '-')
        df = df[f1 & f2 & f3 & f4 & f5]
        return df

    def create_maf(self):

        for center in CENTERS:

            df = pd.DataFrame(columns=maf_columns)
            genomic_df = self.genie_dict[center]['init_genomic_df']

            print '## Input contains %d rows\n' % len(genomic_df.index)
            df.Tumor_Sample_Barcode = genomic_df.Tumor_Sample_Barcode
            df.Center = center.upper()
            df.NCBI_Build = 'GRCh37'
            df.Chromosome = genomic_df.Chromosome
            df.Start_Position = genomic_df.Start_Position
            df.End_Position = genomic_df.End_Position
            df.Strand = genomic_df.Strand.apply(lambda x: '+' if pd.isnull(x) else x)
            df.Reference_Allele = genomic_df.Reference_Allele
            df.Tumor_Seq_Allele1 = genomic_df.Tumor_Seq_Allele1
            df.Tumor_Seq_Allele2 = genomic_df.Tumor_Seq_Allele2
            df = df.fillna(value='', axis=1)

            df = self.remove_indels(df)
            self.write_results(df, center)


def main(args):

    run = Genomic(args)
    run.load_data()
    run.remove_non_breast_cancer_samples()
    run.create_maf()

if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-maf', dest='in_genomic', required=True,
                        help='Path to input genomic genomics file from GENIE')

    parser.add_argument('--in-clinical', dest='in_clinical', required=True,
                        help='The output directory specified by prep_genie_clincial.py is needed in order to'
                             'remove non breast cancer samples from the GENIE genomic file.')

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
