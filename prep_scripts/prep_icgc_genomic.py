#!/usr/bin/env python

"""
SYNOPSIS
     Prepare the ICGC BRCA-UK simple_somatic_mutation.open.BRCA.UK.tsv file into a .maf format
     that is suitable for input into the Hotspots 'Mutation_Counts' pipeline.

NOTES
    - Data was downloaded from https://dcc.icgc.org/releases/current/Projects/BRCA-UK.

EXAMPLES
    ```
    ./prep_icgc_genomic.py \
        --in-vcf ./simple_somatic_mutation.open.BRCA.UK.tsv \
        --out-maf ./icgc_brca_UK.r1.vep.maf
    ```


AUTHOR
    Zachary Zwiesler <zwiesler@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import json
import argparse

import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

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


class VcfToMaf:
    """Contains main steps involved with parsing ICGC data to .maf"""

    def __init__(self, args):
        self.args = args

        self.maf_df = None
        self.icgc_df = None
        self.gene_name_dict = None

        self.num_rows = None
        self.num_samples = None

        self.strand_dict = {'1': '+', '0': '-'}

        # iterator used for logging
        self.i = 0

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'
        self.icgc_df = pd.read_csv(self.args.in_genomics, sep='\t')
        self.maf_df = pd.DataFrame(columns=MAF_columns, index=self.icgc_df.index)

        # these stats are appended to the head of the output .maf file
        self.num_rows = len(self.icgc_df.index)
        self.num_samples = len(self.icgc_df.icgc_sample_id.unique())

    def write_results(self):
        """Write maf dataframe to file."""

        # .maf file
        with open(self.args.out_maf, 'w') as ff:
            ff.write('#version 1.0\n')
            ff.write('#Original .vcf contained %d rows and %d samples\n' % (self.num_rows, self.num_samples))
            self.maf_df.to_csv(ff, sep='\t', index=False)

    def convert_strand(self, strand):
        """Map strand from 1/0 to +/-"""
        return self.strand_dict[str(strand)]

    def create_maf(self):
        """
        Parse each of the following columns from the .vcf file to the .maf file:

        gene_affected --> Hugo_Symbol
        gene_affected --> Entrez_Gene_Id
        project_code --> Center
        assembly_version --> NCBI_Build
        chromosome --> Chromosome
        chromosome_start --> Start_Position
        chromosome_end --> End_Position
        chromosome_strand --> Strand
        reference_genome_allele -> Reference_Allele
        mutated_from_allele -> Tumor_Seq_Allele1
        mutated_to-allele -> Tumor_Seq_Allele2
        icgc_sample_id -> Tumor_Sample_Barcode
        matched_icgc_sample_id -> Matched_Norm_Sample_Barcode
        """

        print '## Input contains %d rows\n' % len(self.icgc_df.index)
        self.maf_df.Center = self.icgc_df.project_code
        self.maf_df.NCBI_Build = self.icgc_df.assembly_version
        self.maf_df.Chromosome = self.icgc_df.chromosome
        self.maf_df.Start_Position = self.icgc_df.chromosome_start
        self.maf_df.End_Position = self.icgc_df.chromosome_end
        self.maf_df.Strand = self.icgc_df.chromosome_strand.apply(lambda x: self.convert_strand(x))
        self.maf_df.Reference_Allele = self.icgc_df.reference_genome_allele
        self.maf_df.Tumor_Seq_Allele1 = self.icgc_df.mutated_from_allele
        self.maf_df.Tumor_Seq_Allele2 = self.icgc_df.mutated_to_allele
        self.maf_df.Tumor_Sample_Barcode = self.icgc_df.icgc_sample_id

    def remove_indels(self):
        """
        As of 05-19-2017, the decision was made to drop all indels from downstream analyses.
        This function does this by removing all rows where the lengths of the Reference_Allele,
        Tumor_Seq_Allele1, and Tumor_Seq_Allele2 are not equal. Note that this allows for any number
        of nucleotide polymorphisms, not just SNPs.
        """

        f1 = (len(self.maf_df.Reference_Allele) == len(self.maf_df.Tumor_Seq_Allele1))
        f2 = (len(self.maf_df.Reference_Allele) == len(self.maf_df.Tumor_Seq_Allele2))
        f3 = (self.maf_df.Reference_Allele != '-')
        self.maf_df = self.maf_df[f1 & f2 & f3]
        print '\n## Output contains %d rows' % len(self.maf_df.index)


def main(args):

    run = VcfToMaf(args)
    run.load_data()
    run.create_maf()
    run.remove_indels()
    run.write_results()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-vcf', dest='in_genomics', required=True,
                        help='Path to input ICGC Genomics Data.')

    parser.add_argument('--out-maf', dest='out_maf', required=True,
                        help='Path to the output MAF file.')
    args = parser.parse_args()

    print '\n## {0}\n## Specified Input:\n{1}\n## {0}'.format(
        '-' * 50, json.dumps(vars(args), indent=4))

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
