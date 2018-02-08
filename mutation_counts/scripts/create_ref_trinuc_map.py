#!/usr/bin/env python2

"""
SYNOPSIS

    Creates a mapping of each transcript-codon-frame key to the corresponding
    reference trinucleotide for all transcripts found in the MAF file.

WARNING

    This script contains PROTOTYPE CODE. Need to significantly test, refactor,
    clean-up and document this code before Publication.

EXAMPLES

    # Quick Test
    ./create_ref_trinuc_map.py \
        --in_maf ../examples/inputs/public/metabric/metabric.subset_genes.vep.maf \
        --gene_annotation_gtf ../pipeline_resources/gene_annotation/r1/subset_genes.CDS_only.Homo_sapiens.GRCh37.85.r1.gtf \
        --ref_fasta ../pipeline_resources/ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        --out_map private/extract_sequence/test_runs/codon2trinuc_map_resultsminimal.codon2trinuc_map.txt

    # DFCI + METABRIC
    ./create_ref_trinuc_map.py \
        --in_maf private/dfci_plus_metabric.r6.filtered_clinical.vep.maf \
        --gene_annotation_gtf ../pipeline_resources/gene_annotation/r1/CDS_only.Homo_sapiens.GRCh37.85.r1.gtf \
        --ref_fasta ../pipeline_resources/ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        --out_map private/extract_sequence/test_runs/dfci_plus_metabric.r6.filtered_clinical.codon2trinuc_map.txt

    # DFCI + Public
    ./create_ref_trinuc_map.py \
        --in_maf private/merged.dfci_plus_public.r6.filtered_clinical.vep.maf \
        --gene_annotation_gtf ../pipeline_resources/gene_annotation/r1/CDS_only.Homo_sapiens.GRCh37.85.r1.gtf \
        --ref_fasta ../pipeline_resources/ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        --out_map private/extract_sequence/test_runs/merged.dfci_plus_public.r6.filtered_clinical.codon2trinuc_map.txt

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
"""

import sys
import os
import time
import copy
import argparse
from pyfaidx import Fasta

import pandas as pd

pd.set_option('display.precision', 2)
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 2000)

VERSION = 'mc_v0.9_dev'

def main(options):

    maf_df = import_annotated_maf(options.in_maf)

    ref_fasta_fai_file = options.ref_fasta + '.fai'

    trinucleotides_map = run_wrapper(maf_df,
                                     options.gene_annotation_gtf,
                                     options.ref_fasta,
                                     ref_fasta_fai_file,
                                     VERSION)

    write_output_table(trinucleotides_map, options.out_map)

def write_output_table(trinucleotides_map, outfile):
    """Write the data in the trinucleotides_map to the specified outfile."""

    fout = open(outfile , 'w')
    header_columns = ["Gene_Symbol", "Transcript_ID", "Codon_Pos", "Strand",
                      "Frame", "Chrom", "Position", "Trinucleotide"]
    fout.write("\t".join(header_columns) + '\n')


    out_header_cols = ['Gene_Symbol', 'Transcript_ID', 'Codon_Pos',
                       'Num_Samples_with_Mutation', 'PValue', 'QValue',
                       'LogLikelihoodRatio', 'LogLikelihood_Alt', 'LogLikelihood_Null']

    for transcript_id in sorted(trinucleotides_map):
        for codon_pos in sorted(trinucleotides_map[transcript_id]):
            frames = sorted(trinucleotides_map[transcript_id][codon_pos])
            assert frames == [0,1,2]
            for frame in frames:

                out_data = trinucleotides_map[transcript_id][codon_pos][frame]
                cols = list()
                cols.append(out_data['Gene_Symbol'])
                cols.append(transcript_id)
                cols.append(str(codon_pos))
                cols.append(out_data['Strand'])
                cols.append(str(frame))
                cols.append(out_data['Chrom'])
                cols.append(str(out_data['Position']))
                cols.append(out_data['Trinucleotide'])
                fout.write("\t".join(cols) + "\n")

    fout.close()

def run_wrapper(maf_df, gene_annotation_gtf_file, ref_fasta_file, ref_fasta_fai_file, version_check):
    """Wrapper function to create a mapping of each transcript-codon-frame key
    to the corresponding reference trinucleotide for all transcripts found in the
    input MAF file."""

    print "##"
    print "## Enter %s's run_wrapper()." % os.path.basename(__file__)
    print "##"

    # Consistency checks
    if ref_fasta_file + '.fai' != ref_fasta_fai_file:
        raise Exception("Unexpected ref_fasta_fai location! Please put the file in same directory as the ref_fasta_file.")

    if version_check != VERSION:
        # This check is to ensure that all deployed scripts inside IntelCCC;s Reliance Point are pulled from the same
        # version of the Git repo. Very unelegant solution but necessary due to limitation of Reliance Point.
        print "##"
        print "## ERROR: Detected Version Incompatibility inside %s's run_wrapper function." % os.path.basename(__file__)
        print "## ERROR: Specified version_check:", repr(version_check)
        print "## ERROR: Script's VERSION:", repr(VERSION)
        raise Exception("Detected Version Incompatibility inside %s's run_wrapper function." % os.path.basename(__file__))


    transcript_ids = sorted((maf_df['Transcript_ID'].unique()))

    transcript_id2gtf_data_list, transcript_id2gene = import_gene_annotation_gtf(gene_annotation_gtf_file)

    ref_genome = Fasta(ref_fasta_file)

    ref_trinucleotides_map = create_ref_trinucleotides_map(transcript_ids,
                                                           transcript_id2gtf_data_list,
                                                           transcript_id2gene,
                                                           ref_genome)

    print "##"
    print "## Exit %s's run_wrapper()." % os.path.basename(__file__)

    return ref_trinucleotides_map

def create_ref_trinucleotides_map(transcript_ids,
                                      transcript_id2gtf_data_list,
                                      transcript_id2gene,
                                      ref_genome):
    """Create a mapping from each transcript-codon-frame key to the corresponding
    reference trinucleotide."""

    ref_trinucleotides_map = dict()

    success_transcripts = list()
    missing_transcripts = list()
    problem_transcripts = list()

    for transcript_id in transcript_ids:

        transcript_id = str(transcript_id)

        if transcript_id == 'nan':
            continue

        if not transcript_id.startswith('ENST'):
            continue

        if transcript_id not in transcript_id2gtf_data_list:
            print "## WARNING: Encounter missing transcript (%s)" % transcript_id
            missing_transcripts.append(transcript_id)
            continue

        per_transcript_gtf_data_list = transcript_id2gtf_data_list[transcript_id]

        codonpos2trinucleotides = map_codons2trinucleotides(transcript_id2gtf_data_list,
                                                            transcript_id, ref_genome)

        if codonpos2trinucleotides == None:
            problem_transcripts.append(transcript_id)
            continue


        success_transcripts.append(transcript_id)

        assert transcript_id not in ref_trinucleotides_map
        ref_trinucleotides_map[transcript_id] = dict()

        for codon_pos in sorted(codonpos2trinucleotides):

            assert codon_pos not in ref_trinucleotides_map[transcript_id]
            ref_trinucleotides_map[transcript_id][codon_pos] = dict()

            frames = sorted(codonpos2trinucleotides[codon_pos])

            assert frames  == [0,1,2]

            for frame in frames:
                gene_symbol = transcript_id2gene[transcript_id]
                strand = codonpos2trinucleotides[codon_pos][frame]['strand']
                chrom = codonpos2trinucleotides[codon_pos][frame]['chrom']
                position = codonpos2trinucleotides[codon_pos][frame]['position'] + 1  # Convert from 0-base to 1-base
                trinucleotide = codonpos2trinucleotides[codon_pos][frame]['trinucleotide']

                assert frame not in ref_trinucleotides_map[transcript_id][codon_pos]
                ref_trinucleotides_map[transcript_id][codon_pos][frame] = dict()
                ref_trinucleotides_map[transcript_id][codon_pos][frame]['Gene_Symbol'] = gene_symbol
                ref_trinucleotides_map[transcript_id][codon_pos][frame]['Strand'] = strand
                ref_trinucleotides_map[transcript_id][codon_pos][frame]['Chrom'] = chrom
                ref_trinucleotides_map[transcript_id][codon_pos][frame]['Position'] = position
                ref_trinucleotides_map[transcript_id][codon_pos][frame]['Trinucleotide'] = trinucleotide

    print "##"
    print "## Successfully processed %s transcripts." % len(success_transcripts)

    if len(missing_transcripts) > 0:
        print "##"
        print "## DEBUG: Encountered %d missing transcripts." % len(missing_transcripts)

    if len(problem_transcripts) > 0:
        print "##"
        print "## DEBUG: Encountered and skipped %d problem transcripts." % len(problem_transcripts)

    return ref_trinucleotides_map

def invert_codonpos2trinucleotides(codonpos2trinucleotides):
    """Invert the codonpos2trinucleotides to account for the fact that the
    transcript is located on the negative (-) strand."""

    protein_length = max(codonpos2trinucleotides)

    inverted_codonpos2trinucleotides = dict()

    for codonpos in codonpos2trinucleotides:
        inverted_codonpos = (protein_length + 1) - codonpos
        assert inverted_codonpos not in inverted_codonpos2trinucleotides
        inverted_codonpos2trinucleotides[inverted_codonpos] = dict()

        for frame in codonpos2trinucleotides[codonpos]:

            assert frame in [0,1,2]
            inverted_frame = 2 - frame

            assert inverted_frame not in inverted_codonpos2trinucleotides[inverted_codonpos]

            inverted_codonpos2trinucleotides[inverted_codonpos][inverted_frame] = codonpos2trinucleotides[codonpos][frame]

    return inverted_codonpos2trinucleotides

def map_codons2trinucleotides(transcript_id2gtf_data_list, transcript_id, ref_genome):
    """Map each codon in the specified transcript to corresponding trinucleotides
    (around the 1st, 2nd, and 3rd base of the codon)

    Note
    ----
    (1) This function will first number the codons assuming that transcript coding
        direction is in the forward strand. If the strand is reverse strand,
        the codon_pos is then inverted afterwards.
    """

    codonpos2trinucleotides = dict()

    gtf_data_list = transcript_id2gtf_data_list[transcript_id]

    strand = set([x['strand'] for x in gtf_data_list])
    assert len(strand) == 1
    strand = list(strand)[0]
    assert strand in ['-', '+']

    # Sort
    if strand == '+':
        gtf_data_list = sorted(gtf_data_list, key=lambda k: k['exon_number'])
    else:
        gtf_data_list = sorted(gtf_data_list, key=lambda k: k['exon_number'], reverse=True)

    codon_pos = 1
    curr_frame = 0

    # For consistency check
    transcript_cds_seq = ""

    # TODO: Check that the exons in gtf_data_list are contiguous.
    # TODO: Check that 'chrom' of all gtf_data are identical

    for gtf_data in gtf_data_list:

        # Consistency checks
        assert gtf_data['feature'] == 'CDS'
        assert gtf_data['transcript_id'] == transcript_id

        chrom = gtf_data['seqname']
        start_coord = gtf_data['start'] - 1  # Convert to 0-base
        end_coord = gtf_data['end'] - 1  # Convert to 0-base

        if start_coord > end_coord:
            print "start_coord:", start_coord
            print "end_coord:", end_coord
            raise Exception("start_coord > end_coord")


        transcript_cds_seq += str(ref_genome[chrom][start_coord:end_coord+1])

        index = 0

        # Temporary dictionary containing mappings from coordinate to nucleotide
        # letter ('A', 'C', 'T', 'G')
        offset = 2  # Handle codons near exon boundaries.
        coordinate2nucleotide = dict()
        dna_sequence_with_offset = str(ref_genome[chrom][start_coord-offset:end_coord+1+offset])
        for coord in xrange(start_coord-offset, end_coord+offset):
            assert coord not in coordinate2nucleotide
            nucleotide_letter = dna_sequence_with_offset[index]
            assert nucleotide_letter in ['A', 'C', 'T', 'G']
            coordinate2nucleotide[coord] = nucleotide_letter
            index += 1

        for coord in range(start_coord, end_coord+1):

            trinucleotide = (coordinate2nucleotide[coord-1] +
                             coordinate2nucleotide[coord] +
                             coordinate2nucleotide[coord+1])

            trinucleotide = normalize_trinucleotide(trinucleotide)

            if codon_pos not in codonpos2trinucleotides:
                codonpos2trinucleotides[codon_pos] = dict()

            codonpos2trinucleotides[codon_pos][curr_frame] = dict()
            codonpos2trinucleotides[codon_pos][curr_frame]['strand'] = strand
            codonpos2trinucleotides[codon_pos][curr_frame]['chrom'] = chrom
            codonpos2trinucleotides[codon_pos][curr_frame]['position'] = coord
            codonpos2trinucleotides[codon_pos][curr_frame]['trinucleotide'] = trinucleotide

            curr_frame += 1
            if curr_frame == 3:
                # Increment to next codon
                curr_frame = 0
                codon_pos += 1

    # Consistency checks that everything is in-frame
    if len(transcript_cds_seq) % 3 != 0:
        print "## WARNING: Encounter Problem Transcript_ID (%s)." % transcript_id,
        print "Problem Transcript CDS Length (%s) is not divisible by 3." % len(transcript_cds_seq)
        return None

    assert curr_frame == 0

    protein_length = int(len(transcript_cds_seq) / 3)
    assert (codon_pos - 1) == protein_length

    if strand == "-":
        codonpos2trinucleotides = invert_codonpos2trinucleotides(codonpos2trinucleotides)

    return codonpos2trinucleotides

def normalize_trinucleotide(trinucleotide):
    """Choose between the current or reverse complement representation of
    the trinucleotide (whichever has the lower lexicographical order).
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

def import_gene_annotation_gtf(gene_annotation_gtf):
    """Import the data in the Gene Annotation GTF file into a dictionary object
    with the 'Transcript_ID' as its key."""

    transcript_id2gtf_data_list = dict()

    transcript_id2gene = dict()

    line_num = 0

    for line in open(gene_annotation_gtf, 'rU'):

        line_num += 1

        if line_num % 100000 == 0:
            print "## PROGRESS: Processing GTF line_num: %s." % line_num

        if line.startswith("#"):
            continue

        cols = line.rstrip('\n').split('\t')

        assert len(cols) == 9

        gtf_data = dict()
        gtf_data["seqname"] = cols[0]
        gtf_data["source"] = cols[1]
        gtf_data["feature"] = cols[2]
        gtf_data["start"] = int(cols[3])
        gtf_data["end"] = int(cols[4])
        gtf_data["score"] = cols[5]
        gtf_data["strand"] = cols[6]
        gtf_data["frame"] = cols[7]

        attributes = cols[8].strip(" ").split(';')

        for attribute in attributes:
            if attribute == "":
                continue

            attribute = attribute.strip(" ")

            assert len(attribute.split(" ")) == 2
            attr_key = attribute.split(" ")[0]
            attr_val = attribute.split(" ")[1]
            assert attr_val.startswith('"') and attr_val.endswith('"')
            attr_val = attr_val.strip('"')

            if attr_key == "exon_number":
                attr_val = int(attr_val)

            # Skip 'tag' key, since this key can appear multiple times.
            if attr_key == 'tag':
                continue

            # Ensure their is no duplicated key
            assert attr_key not in gtf_data

            gtf_data[attr_key] = attr_val

        # Keep only CDS feature
        if gtf_data['feature'] != 'CDS':
            continue

        assert "gene_name" in gtf_data
        assert "exon_number" in gtf_data
        assert "transcript_id" in gtf_data

        key = gtf_data['transcript_id']

        if key not in transcript_id2gtf_data_list:
            transcript_id2gtf_data_list[key] = list()
        transcript_id2gtf_data_list[key].append(gtf_data)

        if key not in transcript_id2gene:
            transcript_id2gene[key] = set()

        transcript_id2gene[key].add(gtf_data["gene_name"])

    # Consistency check
    for transcript_id in transcript_id2gene:
        gene_set = transcript_id2gene[transcript_id]
        assert len(gene_set) == 1
        transcript_id2gene[transcript_id] = list(gene_set)[0]


    return transcript_id2gtf_data_list, transcript_id2gene

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
if __name__ == '__main__':

    print "## Enter %s (%s).\n##" % (os.path.basename(__file__), time.asctime())

    start_time = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument("--in_maf", action="store", required=True,
                        metavar='FILE',
                        help="Path to input MAF file.")

    parser.add_argument("--gene_annotation_gtf", action="store",
                        required=True, metavar='FILE',
                        help="Path to the Gene Annotation GTF resource file.")

    parser.add_argument("--ref_fasta", action="store", required=True, metavar='FILE',
                        help="Path to the Reference Genome Fasta file.")

    parser.add_argument("--out_map", action="store", required=True,
                        metavar='FILE',
                        help="Path to output map file")

    options = parser.parse_args()

    print "##", "-" * 50
    print "## Specified Options:"
    print "##   in_maf:", repr(options.in_maf)
    print "##   gene_annotation_gtf:", repr(options.gene_annotation_gtf)
    print "##   ref_fasta:", repr(options.ref_fasta)
    print "##   out_map:", repr(options.out_map)
    print "##", "-" * 50

    main(options)

    print "##"
    print "## Exit %s" % os.path.basename(__file__),
    print '(%s | total_time = %.3f secs).' % (time.asctime(), time.time() - start_time)

    sys.exit(0)
