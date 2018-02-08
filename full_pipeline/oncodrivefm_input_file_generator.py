#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import sys
import numpy as np
import pandas as pd


def process_header_lines(header_lines):
    h = {}
    for l in header_lines:
        split = l.strip('#').strip("\n").split("=", 1)
        try:
            if not split[0] in h:
                h[split[0]] = [split[1]]
            else:
                h[split[0]] += [split[1]]
        except:
            pass
    return h


def process_info_fields(info_str):
    split = info_str.strip("<>\"").split(': ')[1].split("|")
    info_map = dict(zip(split, list(range(0, len(split)))))        
    return info_map


def extract_records(info_str, canonical=False):
    records = info_str.split(",")
    records = [r.split("|") for r in records]
    if canonical:
        records = [r for r in records if r[24] == 'YES']
    return records


def process_mutation_scores(records, info_fields_map):
    score_map = {
        'stSNV': {
            "SIFT": 0,
            "PPH2": 1,
            "MA": 3.5
        },
        'fsindel': {
            "SIFT": 0,
            "PPH2": 1,
            "MA": 3.5
        },
        'sSNV': {
            "SIFT": 1,
            "PPH2": 0,
            "MA": -2
        }
    }

    consequence_map = {
        "missense_variant": "nsSNV",
        "synonymous_variant": "sSNV",
        "frameshift_variant": "fsindel",
        "stop_gained": "stSNV"
    }

    gene = []
    sift = []
    pph = []
    ma = []
    for r in records:
        if r[info_fields_map['Consequence']] in consequence_map:
            gene.append(r[info_fields_map['Gene']])     
            if r[info_fields_map['Consequence']] == "missense_variant":
                sift.extend(
                    [eval_else(v, float) for v in r[info_fields_map['SIFT_score']].split("&")]
                )
                pph.extend(
                    [eval_else(v, float) for v in r[info_fields_map['Polyphen2_HDIV_score']].split("&")]
                )
                ma.extend(
                    [eval_else(v, float) for v in r[info_fields_map['MutationAssessor_score']].split("&")]
                )
            else:
                sift.append(
                    score_map[consequence_map[r[info_fields_map['Consequence']]]]['SIFT']
                )
                pph.append(
                    score_map[consequence_map[r[info_fields_map['Consequence']]]]['PPH2']
                )
                ma.append(
                    score_map[consequence_map[r[info_fields_map['Consequence']]]]['MA']
                )
        else:
            break

    if len(gene) >= 1:
        if all(x == gene[0] for x in gene):
            gene = gene[0]
        else:
            print("[WARNING] multiple gene ids are associated with this variant:", records, file=sys.stderr)
            gene = max(set(gene), key=gene.count)
    else:
        gene = None

    return [gene, eval_else(sift, min), eval_else(pph, max), eval_else(ma, max)]


def eval_else(v, func, ret=np.nan):
    try:
        return func(v)
    except:
        return ret


# vcfFile = "/Users/strucka/Projects/dockerized_tools/vep/cromwell-executions/test/6606e7a3-39dc-4556-a86a-5a603e777cd0/call-variant_effect_predictor/vep_output.vcf"
def main(vcf_file, output_file):
    # read files
    header_lines = os.popen('head -5000 ' + vcf_file).readlines()
    header_lines = [l for l in header_lines if l.startswith('#')]
    vcf_header = process_header_lines(header_lines)
    info_fields = process_info_fields(vcf_header['INFO'][0])

    vcf = pd.read_table(vcf_file, header=len(header_lines)-1, na_values="./.:.:.", low_memory=False)

    # extract sift, pph2 and ma scores
    info_df = vcf['INFO'].apply(lambda x: extract_records(x, canonical=False)).apply(lambda x: process_mutation_scores(x, info_fields)).apply(pd.Series)
    info_df.columns = ['Gene', 'SIFT_score', 'Polyphen2_HDIV_score', 'MutationAssessor_score']

    # merge with vcf
    min_vcf = vcf.drop(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"], axis=1, inplace=False)
    expanded_vcf = pd.concat(
        [
            min_vcf.drop(["INFO"], axis=1, inplace=False),
            info_df
        ],
        axis=1
    )

    # reshape
    melted = pd.melt(expanded_vcf,
                     id_vars=['Gene',
                              'SIFT_score',
                              'Polyphen2_HDIV_score',
                              'MutationAssessor_score'],
                     var_name="sample",
                     value_name="mutation_status")

    # drop na, reorder columns and rename columns
    filtered = melted.dropna()[melted['sample'] != "NORMAL"]
    reordered = filtered[["sample", "Gene", "SIFT_score", "Polyphen2_HDIV_score", "MutationAssessor_score"]]
    reordered.columns = ["SAMPLE", "GENE", "SIFT", "PPH2", "MA"]

    # write to output file
    reordered.to_csv(output_file, sep="\t", header=True, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("vcf_file",                        
                        type=str,
                        help="VEP annotated VCF file with cols: Gene, SIFT_score, Polyphen2_HDIV_score, MutationAssessor_score")
    parser.add_argument("-o",
                        dest="output_file",
                        type=str,
                        default="oncodrivefm_input.tdm",
                        help="output file")
    args = parser.parse_args()

    main(args.vcf_file, args.output_file)
