task merge {
    Array[File]+ vcf_files
    Array[File]+ vcf_files_tbi
    File? use_header
    File? regions_file
    Array[String]? apply_filters
    Array[String]? info_rules
    Array[String]? regions
    Boolean force_samples = true
    Boolean print_header = false
    Boolean no_version = false
    String? merge
    String output_type = 'v'
    String output_file_prefix
    Int threads = 3

    command {
        # Handling the optional arrays
        if [ -n '${sep="," apply_filters}' ]; then
            FILTER_FLAG=--apply-filters ${sep="," apply_filters}
        else
            FILTER_FLAG=''
        fi

        if [ -n '${sep="," info_rules}' ]; then
            RULES_FLAG=--info-rules ${sep="," info_rules}
        else
            RULES_FLAG=''
        fi

        if [ -n '${sep="," regions}' ]; then
            REGIONS_FLAG=--regions ${sep="," regions}
        else
            REGIONS_FLAG=''
        fi

        bcftools merge \
        ${true="--force-samples" false="" force_samples} \
        ${true="--print-header" false="" print_header} \
        ${true="--no-version" false="" no_version} \
        $FILTER_FLAG \
        $RULES_FLAG \
        $REGIONS_FLAG \
        ${"--regions-file " + regions_file} \
        ${"--use-header " + use_header} \
        ${"--merge " + merge} \
        ${"--threads " + threads} \
        --output-type ${output_type} \
        --output ${output_file_prefix}.vcf \
        ${sep=" " vcf_files}
    }

    output {
        File merged_vcf = "${output_file_prefix}.vcf"
    }

    runtime {
        docker: "bcftools"
    }

}

task chromNamesToEnsembl {
    File inputFile
    String ext = "vcf"
    String outputFilePrefix

    command {
        sed 's/^chr//g' ${inputFile} > ${outputFilePrefix}_ensemblChroms.${ext}
    }

    output {
        File file_with_ensemblchroms = "${outputFilePrefix}_ensemblchroms.${ext}"
    }

    runtime {
        docker: "vcf2maf"
    }
}

task variant_effect_predictor {
    File inputFile
    String outputFileName
    File refFasta
    File refFastaFai
    File cacheDir
    File? pluginsDir
    Int? fork

    String species = "homo_sapiens"
    String assembly = "GRCh37"
    String format = "vcf"

    # shortcut for common flags
    Boolean? everything

    # output options
    Boolean? humdiv
    Boolean? gene_phenotype
    Boolean? regulatory
    Boolean? phased
    Boolean? allele_number
    Boolean? total_length
    Boolean? numbers
    Boolean? domains
    Boolean? no_escape
    Boolean? keep_csq
    Boolean? no_consequences
    Boolean? variant_class

    String? cell_type
    String? individual
    String? sift
    String? polyphen
    String? vcf_info_field
    String? terms

    # identifiers
    Int? shift_hgvs
    Boolean? hgvs
    Boolean? protein
    Boolean? symbol
    Boolean? ccds
    Boolean? uniprot
    Boolean? tsl
    Boolean? appris
    Boolean? canonical
    Boolean? biotype
    Boolean? xref_refseq

    # co-located variants
    Boolean? check_existing
    Boolean? check_alleles
    Boolean? check_svs
    Boolean? gmaf
    Boolean? maf_1kg
    Boolean? maf_esp
    Boolean? maf_exac
    Boolean? old_maf
    Boolean? pubmed
    Int? failed

    # output format options
    Boolean? vcf
    Boolean? tab
    Boolean? json
    Boolean? gvcf
    Boolean? minimal

    command {
        dbNSFP_PLUGIN="--plugin dbNSFP,${cacheDir}/dbNSFP.gz,SIFT_score,Polyphen2_HDIV_score,MutationAssessor_score"
        ExAC_PLUGIN="--plugin ExAC,${cacheDir}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"

        variant_effect_predictor.pl \
        --no_progress \
        --no_stats \
        --offline \
        --input_file ${inputFile} \
        --format ${format} \
        --species ${species} \
        --assembly ${assembly} \
        --fasta ${refFasta} \
        --dir_cache ${cacheDir} \
        ${"--dir_plugins " + pluginsDir} \
        ${"--fork " + fork} \
        ${true="--everything" false="" everything} \
        $dbNSFP_PLUGIN \
        $ExAC_PLUGIN \
        ${"--cell_type " + cell_type} \
        ${"--individual " + individual} \
        ${"--sift " + sift} \
        ${"--polyphen " + polyphen} \
        ${"--vcf_info_field " + vcf_info_field} \
        ${"--terms " + terms} \
        ${true="--humdiv" false="" humdiv} \
        ${true="--gene_phenotype" false="" gene_phenotype} \
        ${true="--regulatory" false="" regulatory} \
        ${true="--phased" false="" phased} \
        ${true="--allele_number" false="" allele_number} \
        ${true="--total_length" false="" total_length} \
        ${true="--numbers" false="" numbers} \
        ${true="--domains" false="" domains} \
        ${true="--no_escape" false="" no_escape} \
        ${true="--keep_csq" false="" keep_csq} \
        ${true="--no_consequences" false="" no_consequences} \
        ${true="--variant_class" false="" variant_class} \
        ${"--shift_hgvs " + shift_hgvs} \
        ${true="--hgvs" false="" hgvs} \
        ${true="--protein" false="" protein} \
        ${true="--symbol" false="" symbol} \
        ${true="--ccds" false="" ccds} \
        ${true="--uniprot" false="" uniprot} \
        ${true="--tsl" false="" tsl} \
        ${true="--appris" false="" appris} \
        ${true="--canonical" false="" canonical} \
        ${true="--biotype" false="" biotype} \
        ${true="--xref_refseq" false="" xref_refseq} \
        ${true="--check_existing" false="" check_existing} \
        ${true="--check_alleles" false="" check_alleles} \
        ${true="--check_svs" false="" check_svs} \
        ${true="--gmaf" false="" gmaf} \
        ${true="--maf_1kg" false="" maf_1kg} \
        ${true="--maf_esp" false="" maf_esp} \
        ${true="--maf_exac" false="" maf_exac} \
        ${true="--old_maf" false="" old_maf} \
        ${true="--pubmed" false="" pubmed} \
        ${"--failed " + failed} \
        ${true="--vcf" false="" vcf} \
        ${true="--tab" false="" tab} \
        ${true="--json" false="" json} \
        ${true="--gvcf" false="" gvcf} \
        ${true="--minimal" false="" minimal} \
        --output_file ${outputFileName}
    }

    output {
        File annotatedFile = "${outputFileName}"
    }

    runtime {
        docker: "vep"
    }
}

task vcf2maf {
    File inputVCF
    File? vepAnnotatedInputVCF
    String tmpDir = "."
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai
    File? remapChain
    String ncbiBuild = "GRCh37"
    String species = "homo_sapiens"
    Array[String]? retainInfoCols
    String? tumorId
    String? normalId
    String? vcfTumorId
    String? vcfNormalId
    File? customEnst
    String? mafCenter
    Float? minHomVaf

    String outputFilePrefix

    command {
        if [ -n "${vepAnnotatedInputVCF}" ]; then
            ln -s ${vepAnnotatedInputVCF} ${tmpDir}/$(basename ${inputVCF} .vcf).vep.vcf
        fi

        # retain extra INFO cols
        if [ -n "${sep="," retainInfoCols}" ]; then
           INFOCOLS="--retain-info  ${sep ="," retainInfoCols}"
        else
           INFOCOLS=""
        fi

        perl /home/vcf2maf.pl --input-vcf ${inputVCF} \
                              --output-maf ${outputFilePrefix}.maf \
                              --vep-data ${vepOfflineCacheDir} \
                              --ref-fasta ${refFasta} \
                              --species ${species} \
                              --ncbi-build ${ncbiBuild} \
                              --tmp-dir ${tmpDir} \
                              ${"--remap-chain " + remapChain} \
                              ${"--maf-center " + mafCenter} \
                              ${"--tumor-id " + tumorId} \
                              ${"--normal-id " + normalId} \
                              ${"--vcf-tumor-id " + vcfTumorId} \
                              ${"--vcf-normal-id " + vcfNormalId} \
                              ${"--custom-enst " + customEnst} \
                              ${"--min-hom-vaf " + minHomVaf} \
                              $INFOCOLS

       rm ${tmpDir}/$(basename ${inputVCF} .vcf).vep.vcf
    }

    output {
        File maf = "${outputFilePrefix}.maf"
    }

    runtime {
        docker: "vcf2maf"
    }
}

task taylor {
     File inputMAF
     File refFasta
     File refFastaFai
     File? geneQuery
     File? align100mer
     File? align24mer
     String homopolymer = "TRUE"
     String filterCenterBias = "FALSE"
     String outputFilePrefix = "taylor"

     command <<<
         # link ref files to expected location
         ln -s ${refFasta} /mnt/Homo_sapiens.GRCh37.dna.primary_assembly.fa;
         ln -s ${refFastaFai} /mnt/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai;

         # replace NA strings with what the algorithm expects and remove chr
         awk -F"\t" -v OFS="\t" '{for (i=1; i<=NF; i++) { if ($i == ".") $i="NA"} print $0}' ${inputMAF} | sed 's#chr##' > tmp.maf;

         # call hotspot algorithm
         Rscript /home/hotspot_algo.R --input-maf=tmp.maf \
                                --rdata=/home/hotspot_algo.Rdata \
                                ${"--gene-query=" + geneQuery} \
                                ${"--align100mer=" + align100mer} \
                                ${"--align24mer=" + align24mer} \
                                --homopolymer=${homopolymer} \
                                --filter-centerbias=${filterCenterBias} \
                                --output-file=${outputFilePrefix}.txt;

        # cleanup
        rm tmp.maf
     >>>

     output {
         File hotspots = "${outputFilePrefix}.txt"
     }

     runtime {
         docker: "taylor-lab-hotspots"
     }
}

task oncodrivefm_input_gen {
     File vcfFile
     String outputFilePrefix = "oncodrive_inputgen"

     command <<<
         python <<CODE
         from __future__ import print_function
         
         import os
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
                     print("[WARNING] multiple gene ids are associated with this variant:", r )
                     gene = max(set(gene), key=gene.count)
             else:
                 gene = None
         
             return [gene, eval_else(sift, min), eval_else(pph, max), eval_else(ma, max)]
         
         
         def eval_else(v, func, ret=np.nan):
             try:
                 return func(v)
             except:
                 return ret
         
         
         vcfFile = "${vcfFile}"
         
         # read files
         header_lines = os.popen('head -5000 ' + vcfFile).readlines()
         header_lines = [l for l in header_lines if l.startswith('#')]
         vcf_header = process_header_lines(header_lines)
         info_fields = process_info_fields(vcf_header['INFO'][0])
         
         vcf = pd.read_table(vcfFile, header=len(header_lines)-1, na_values="./.:.:.", low_memory=False)
         
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
         reordered.to_csv("${outputFilePrefix}.tdm", sep="\t", header=True, index=False)
         CODE
     >>>

     output {
         File oncodrivefm_input = "${outputFilePrefix}.tdm"
     }

     runtime {
         docker: "python-numpy-pandas-scikit-learn"
     }
}

task oncodrivefm {
    File inputData
    File pathwayMappingFile
    String estimator
    Int? numSamplings
    Int? pathwayThreshold
    Int? geneThreshold
    Int? cores
    String? outputFormat
    String outputFilePrefix = "oncodrive"

    command {
        oncodrivefm \
        -e ${estimator} \
        -m ${pathwayMappingFile} \
        ${"-N " + numSamplings} \
        ${"--gt " + geneThreshold} \
        ${"--pt " + pathwayThreshold} \
        ${"-j " + cores} \
        ${"--output-format " + outputFormat} \
        -n ${outputFilePrefix} \
        ${inputData}
    }

    output {
        File genes = "${outputFilePrefix}-genes.tsv"
        File pathways = "${outputFilePrefix}-pathways.tsv"
    }

    runtime {
        docker: "oncodrivefm"
    }
}


workflow hotspot {
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai

    # merge the VCF files
    call merge

    call chromNamesToEnsembl {
        input:  inputFile = merge.merged_vcf
    }

    call variant_effect_predictor {
        input: inputFile = chromNamesToEnsembl.file_with_ensemblchroms,
               cacheDir = vepOfflineCacheDir,
               refFasta = refFasta,
               refFastaFai = refFastaFai
    }

    ## ONCODRIVE ##
    call oncodrivefm_input_gen {
        input: vcfFile = variant_effect_predictor.annotatedFile
    }

    call oncodrivefm {
        input: inputData = oncodrivefm_input_gen.oncodrivefm_input
    }

    ## Taylor ##
    call vcf2maf {
        input: inputVCF = variant_effect_predictor.annotatedFile,
               vepAnnotatedInputVCF = variant_effect_predictor.annotatedFile,
               vepOfflineCacheDir = vepOfflineCacheDir,
               refFasta = refFasta,
               refFastaFai = refFastaFai
    }

    call taylor {
        input: inputMAF = vcf2maf.maf,
               refFasta = refFasta,
               refFastaFai = refFastaFai
    }
}
