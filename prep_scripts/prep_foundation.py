#!/usr/bin/env python2

"""
SYNOPSIS
    Prepares the Foundation data into a format that is suitable for input
    into the Hotspots 'Mutation_Counts' pipeline.

EXAMPLES

    ./prep_foundation.py \
        --in-clinical ./foundation_one_breast_tumors.txt \
        --out-clinical ./foundation.clinical.r1.tsv \
        --out-genomic ./foundation.r1.maf

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

OUT_COLS = [
    'Tumor_Sample_Barcode',
    'Center',
    'ER_Status',
    'PR_Status',
    'HER2_Status',
    'Biopsy_Site_Type',
    'Histology_Type',
    'Gender'
]

MAF_columns = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
               'Chromosome', 'Start_Position', 'End_Position', 'Strand',
               'Variant_Classification', 'Variant_Type', 'Reference_Allele',
               'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
               'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
               'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2',
               'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',
               'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2',
               'Verification_Status', 'Validation_Status', 'Mutation_Status',
               'Sequencing_Phase', 'Sequence_Source', 'Validation_Method',
               'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID',
               'Matched_Norm_Sample_UUID']

histology_type_dict = {
    'breast invasive ductal carcinoma (idc)': 'Invasive_Ductal_Carcinoma',
    'breast invasive lobular carcinoma (ilc)': 'Invasive_Lobular_Carcinoma',
    'breast metaplastic carcinoma': 'Other_Invasive_Breast_Carcinoma',
    'breast inflammatory carcinoma': 'Other_Invasive_Breast_Carcinoma',
    'breast ductal carcinoma in situ (dcis)': 'Other_Invasive_Breast_Carcinoma',
    'breast carcinoma (nos)': 'Other_Invasive_Breast_Carcinoma',
    'breast adenoid cystic carcinoma': 'Other_Invasive_Breast_Carcinoma',
    'breast assorted': 'Other_Breast_Cancer'
}


class Clinical:

    def __init__(self, args):

        self.args = args

        self.foundation_df = None
        self.clinical_df = pd.DataFrame(columns=OUT_COLS)
        self.maf_df = pd.DataFrame(columns=MAF_columns)

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.foundation_df = pd.read_csv(self.args.in_clinical, sep='\t')

    def write_results(self):
        """Write clinical dataframe to file."""

        # clinical file
        with open(self.args.out_clinical, 'w') as ff:
            ff.write('#version 1.0\n')
            self.clinical_df.to_csv(ff, sep='\t', index=False)

        # genomic file
        with open(self.args.out_genomic, 'w') as ff:
            ff.write('#version 1.0\n')
            self.maf_df.to_csv(ff, sep='\t', index=False)

    def get_biopsy_site_type(self, sample_id):
        """Return the corresponding biopsy type value that the hotspot pipeline accepts"""

        f1 = (self.foundation_df['sample name'] == sample_id)
        tissue_of_origin = self.foundation_df[f1]['tissue of origin'].tolist()[0]

        if tissue_of_origin == 'Breast':
            return 'Primary'
        else:
            return 'Metastatic'

    def get_histology_type(self, sample_id):
        """Convert oncotree code to histology type"""

        f1 = (self.foundation_df['sample name'] == sample_id)
        disease_ontology = self.foundation_df[f1]['disease ontology'].tolist()[0]

        if disease_ontology in histology_type_dict:
            return histology_type_dict[disease_ontology]
        else:
            print '## WARNING: Sample removed. This disease ontology was not found: %s' % disease_ontology
            return 'removeme'

    @staticmethod
    def parse_allele(allele):

        if pd.isnull(allele) or allele.startswith('cannotParseTranscriptEffect'):
            return 'removeme'
        else:
            return allele.upper()

    @staticmethod
    def fix_end_position_for_non_snps(df):
        if df.Reference_Allele == 'removeme':
            return 0
        return str(int(df.Start_Position) + len(df.Reference_Allele) - 1)

    def create_clinical_tsv(self):

        print '\n## Creating clinical file...'
        print'## Clinical file contains %d samples' % len(self.foundation_df['sample name'].unique())
        self.clinical_df.Tumor_Sample_Barcode = self.foundation_df['sample name'].unique().tolist()
        self.clinical_df.Center = 'FOUNDATION'
        self.clinical_df.ER_Status = 'Unknown'
        self.clinical_df.PR_Status = 'Unknown'
        self.clinical_df.HER2_Status = 'Unknown'
        self.clinical_df.Biopsy_Site_Type = self.clinical_df.Tumor_Sample_Barcode.apply(lambda x: self.get_biopsy_site_type(x))
        self.clinical_df.Histology_Type = self.clinical_df.Tumor_Sample_Barcode.apply(lambda x: self.get_histology_type(x))
        self.clinical_df.Gender = self.foundation_df['patient gender'].apply(lambda x: x.title())
        self.clinical_df = self.clinical_df.fillna(value='Unknown', axis=1)

        # drop non breast cancers
        f1 = (self.clinical_df.Histology_Type != 'removeme')
        self.clinical_df = self.clinical_df[f1]

    def create_genomic(self):

        print '\n## Creating genomic file...'
        print '## Input contains %d rows\n' % len(self.foundation_df['sample name'].unique())
        self.maf_df.Tumor_Sample_Barcode = self.foundation_df['sample name']
        self.maf_df.Center = 'FOUNDATION'
        self.maf_df.NCBI_Build = 'GRCh37'
        self.maf_df.Chromosome = self.foundation_df.chromosome
        self.maf_df.Start_Position = self.foundation_df.position
        self.maf_df.Strand = self.foundation_df.strand.apply(lambda x: '+' if pd.isnull(x) else x)
        self.maf_df.Reference_Allele = self.foundation_df['ref allele'].apply(lambda x: self.parse_allele(x))
        self.maf_df.Tumor_Seq_Allele1 = self.foundation_df['ref allele'].apply(lambda x: self.parse_allele(x))
        self.maf_df.Tumor_Seq_Allele2 = self.foundation_df['alt allele'].apply(lambda x: self.parse_allele(x))
        self.maf_df.End_Position = self.maf_df.apply(lambda x: self.fix_end_position_for_non_snps(x), axis=1)
        self.maf_df = self.maf_df.fillna(value='', axis=1)

        # drop nan alleles
        f1 = (self.maf_df.Reference_Allele != 'removeme')
        f2 = (self.maf_df.Tumor_Seq_Allele1 != 'removeme')
        f3 = (self.maf_df.Tumor_Seq_Allele2 != 'removeme')
        self.maf_df = self.maf_df[f1 & f2 & f3]


def main(args):

    run = Clinical(args)
    run.load_data()
    run.create_clinical_tsv()
    run.create_genomic()
    run.write_results()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-clinical', dest='in_clinical', required=True,
                        help='Path to input Foundation data file')

    parser.add_argument('--out-clinical', dest='out_clinical', required=True,
                        help='Path to output Foundation clinical file')

    parser.add_argument('--out-genomic', dest='out_genomic', required=True,
                        help='Path to output Foundation genomic file')

    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
