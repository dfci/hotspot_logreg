#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the MSK 10k data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

NOTES

   - Public dataset downloaded from cBioPortal datahub:
    https://github.com/cBioPortal/datahub/blob/master/public/msk_impact_2017.tar.gz

EXAMPLES

    ./prep_msk_genomic.py \
        --in-maf ./data_mutations_extended.txt \
        --in-clinical . \
        --out-dir .

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

        self.template_dict = {
            'num_rows': 0,
            'num_samples': 0,
            'keep_these_sample_ids': []
        }
        self.panel_dict = {
            'MSK-IMPACT341': self.template_dict.copy(),
            'MSK-IMPACT410': self.template_dict.copy()
        }

        self.i = 0

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.genomic_df = pd.read_csv(self.args.in_genomic, sep='\t', skiprows=1, low_memory=False)
        for panel in self.panel_dict:
            filename = '%s/msk10k.%s.clinical.r2.tsv' % (self.args.in_clinical, panel.replace('MSK-', '').lower())
            clinical_df = pd.read_csv(filename, sep='\t', skiprows=1)
            self.panel_dict[panel]['num_rows'] = len(clinical_df.index)
            self.panel_dict[panel]['num_samples'] = len(clinical_df.Tumor_Sample_Barcode.unique())
            self.panel_dict[panel]['keep_these_sample_ids'] = clinical_df.Tumor_Sample_Barcode.unique()

    def write_results(self, maf_df, panel):
        """Write clinical dataframe to file."""

        # clinical file
        filename = '%s/msk10k.%s.r4.maf' % (self.args.out_dir, panel.replace('MSK-', '').lower())
        with open(filename, 'w') as ff:
            ff.write('#version 1.0\n')
            ff.write('#Original .maf contained %d rows and %d samples\n' % (
                     self.panel_dict[panel]['num_rows'],
                     self.panel_dict[panel]['num_samples']))
            maf_df.to_csv(ff, sep='\t', index=False)

    def create_maf(self):

        for panel in self.panel_dict:

            print '## Input contains %d rows\n' % self.panel_dict[panel]['num_rows']

            # remove non breast cancer samples
            f1 = (self.genomic_df.Tumor_Sample_Barcode.isin(self.panel_dict[panel]['keep_these_sample_ids']))
            genomic_df = self.genomic_df[f1]

            maf_df = pd.DataFrame(columns=MAF_columns)
            maf_df.Tumor_Sample_Barcode = genomic_df.Tumor_Sample_Barcode
            maf_df.Center = panel
            maf_df.NCBI_Build = 'GRCh37'
            maf_df.Chromosome = genomic_df.Chromosome
            maf_df.Start_Position = genomic_df.Start_Position
            maf_df.End_Position = genomic_df.End_Position
            maf_df.Strand = genomic_df.Strand.apply(lambda x: '+' if pd.isnull(x) else x)
            maf_df.Reference_Allele = genomic_df.Reference_Allele
            maf_df.Tumor_Seq_Allele1 = genomic_df.Tumor_Seq_Allele1
            maf_df.Tumor_Seq_Allele2 = genomic_df.Tumor_Seq_Allele2
            maf_df = maf_df.fillna(value='', axis=1)

            self.write_results(maf_df, panel)


def main(args):

    run = Genomic(args)
    run.load_data()
    run.create_maf()

if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-maf', dest='in_genomic', required=True,
                        help='Path to input genomic .maf file from cBioPortal')

    parser.add_argument('--in-clinical', dest='in_clinical', required=True,
                        help='The output of prep_msk_clincial.py is needed in order to'
                             'remove non breast cancer samples from the MSK genomic file.')

    parser.add_argument('--out-dir', dest='out_dir', required=True,
                        help='Path to output genomic file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
