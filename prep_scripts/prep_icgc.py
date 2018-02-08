#!/usr/bin/env python

"""
SYNOPSIS
    Prepare any ICGC dataset into a MAF format that is suitable for input into
    the Hotspots 'Mutation_Counts' pipeline.

NOTES

    (1) Assume that 'Reference_Allele'/'Tumor_Seq_Allele1'/'Tumor_Seq_Allele2'
        values in input ('--in_genomics') already satisfy the MAF-specification.
    (2) Data was downloaded from https://dcc.icgc.org/releases/current/Projects/BRCA-UK.
        'simple_somatic_mutation.open.BRCA.UK.tsv.gz' was used.

EXAMPLES

    # Extract genomics data into suitable output MAF format.
    ./prep_icgc.py \
        -i ./data/icgc/simple_somatic_mutation.open.BRCA-UK.tsv \
        -o ./data/icgc/icgc_brca_uk.r1.maf \
        --donor-file ./data/icgc/donor.BRCA-UK.tsv \
        --sample-file ./data/icgc/sample.BRCA-UK.tsv \
        --out-clinical ./data/icgc/icgc_brca_uk.clinical.r1.tsv

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

clinical_columns = [
    'Tumor_Sample_Barcode',
    'Center',
    'ER_Status',
    'PR_Status',
    'HER2_Status',
    'Biopsy_Site_Type',
    'Histology_Type',
    'Gender'
]

donor_columns = [
    'icgc_donor_id',
    'project_code',
    'study_donor_involved_in ',
    'submitted_donor_id',
    'donor_sex',
    'donor_vital_status',
    'disease_status_last_followup',
    'donor_relapse_type',
    'donor_age_at_diagnosis',
    'donor_age_at_enrollment',
    'donor_age_at_last_followup',
    'donor_relapse_interval',
    'donor_diagnosis_icd10',
    'donor_tumour_staging_system_at_diagnosis',
    'donor_tumour_stage_at_diagnosis',
    'donor_tumour_stage_at_diagnosis_supplemental',
    'donor_survival_time',
    'donor_interval_of_last_followup',
    'prior_malignancy',
    'cancer_type_prior_malignancy',
    'cancer_history_first_degree_relative'
]

sample_columns = [
    'icgc_sample_id',
    'project_code',
    'submitted_sample_id',
    'icgc_specimen_id',
    'submitted_specimen_id',
    'icgc_donor_id',
    'submitted_donor_id',
    'analyzed_sample_interval',
    'percentage_cellularity',
    'level_of_cellularity',
    'study',
]


class Utilities:
    """
    Contains basic utilities used by all classes.
    This class also contains the data in Pandas dataframes.
    """

    def __init__(self, args):
        self.args = args

        self.maf_df = None
        self.icgc_df = None
        self.donor_df = None
        self.sample_df = None
        self.clinical_df = None
        self.gene_name_dict = None

        self.num_rows = None
        self.num_samples = None
        self.strand_dict = {'1': '+', '0': '-'}

        # iterator used for logging
        self.i = 0

    def convert_strand(self, strand):
        """Map strand from 1/0 to +/-"""
        return self.strand_dict[str(strand)]

    def get_gender(self, donor_id):
        """
        Returns gender from icgc_sample_id

        :param donor_id: ICGC Sample ID
        :return: Male or Female
        """
        f1 = (self.donor_df.icgc_donor_id == donor_id)
        return self.donor_df[f1]['donor_sex'].tolist()[0].title()


class VcfToMaf(Utilities):
    """Contains main steps involved with parsing ICGC data to .maf"""

    def __init__(self, args):
        Utilities.__init__(self, args)

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


class Clinical(Utilities):
    """
    For this project, the only clinical information we have is based off of the
    project name: "Breast Triple Negative/Lobular Cancer". From this,
    we can assume that all samples' hormone receptor status is ER, PR, and HER2
    negative, and histology type is Invasive Lobular Carcinoma.

    Gender was pulled from the donor.${center}.tsv file downloaded from here:
    https://dcc.icgc.org/releases/current/Projects/${center}

    Required columns are:
     - Tumor_Sample_Barcode
     - Center
     - ER_Status
     - PR_Status
     - HER2_Status
     - Biopsy_Site_Type
     - Histology_Type
     - Gender
    """

    def __init__(self, args):
        Utilities.__init__(self, args)

    def load_data(self):
        """Loads all input data from file to Pandas dataframes"""

        print '##\n## Loading data...'

        # Donor file
        if self.args.donor_file:
            if not self.args.sample_file:
                print '## ERROR: To generate the clinical file, you must provide the donor and sample '
                print '##        file paths because the donor ID is mapped to the sample ID in the'
                print '##        mutations file through these.'

            self.clinical_df = pd.DataFrame(columns=clinical_columns)
            self.donor_df = pd.read_csv(self.args.donor_file, sep='\t', header=0)
            self.sample_df = pd.read_csv(self.args.sample_file, sep='\t', header=0)

    def write_results(self):
        """Write clinical dataframe to file."""

        # clinical file
        with open(self.args.out_clinical, 'w') as ff:
            ff.write('#version 1.0\n')
            self.clinical_df.to_csv(ff, sep='\t', index=False)

    def create_clinical_tsv(self):
        """
        Creates clinical file from the three dataframes donor_df, specimen_df, and sample_df,
        as well as assumes hormone receptor status and histology type through the title.
        """

        print '\n## Creating clinical file...'
        print'## Sample file contains %d samples' % len(self.sample_df.index)
        self.clinical_df.Tumor_Sample_Barcode = self.sample_df.icgc_sample_id
        self.clinical_df.Center = self.sample_df.project_code
        self.clinical_df.ER_Status = 'Negative'
        self.clinical_df.PR_Status = 'Negative'
        self.clinical_df.HER2_Status = 'Negative'
        self.clinical_df.Biopsy_Site_Type = 'Primary'
        self.clinical_df.Histology_Type = 'Invasive_Lobular_Carcinoma'
        self.clinical_df.Gender = self.sample_df.icgc_donor_id.apply(lambda x: self.get_gender(x))
        self.clinical_df = self.clinical_df.fillna(value='Unknown', axis=1)



def main(args):

    run = VcfToMaf(args)
    cli = Clinical(args)

    run.load_data()
    cli.load_data()

    run.create_maf()
    run.remove_indels()
    cli.create_clinical_tsv()

    run.write_results()
    cli.write_results()


if __name__ == '__main__':

    print '## Enter %s (%s).\n##' % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='in_genomics', required=True,
                        help='Path to input ICGC Genomics Data.')

    parser.add_argument('-o', dest='out_maf', required=True,
                        help='Path to the output MAF file.')

    parser.add_argument('--donor-file', dest='donor_file', required=False,
                        help='Path to donor clinical file.')

    parser.add_argument('--sample-file', dest='sample_file', required=False,
                        help='Path to sample clinical file.')

    parser.add_argument('--out-clinical', dest='out_clinical', required=False,
                        help='Path to output clinical file.')

    args = parser.parse_args()

    print '##', '-' * 50
    print '## Specified Options:'
    print '##   in_genomics: ', repr(args.in_genomics)
    print '##   out_maf:', repr(args.out_maf)
    print '##   donor_file: ', repr(args.donor_file)
    print '##   sample_file: ', repr(args.sample_file)
    print '##   out_clinical: ', repr(args.out_clinical)
    print '##', '-' * 50

    main(args)

    print '##'
    print '## Exit %s' % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
