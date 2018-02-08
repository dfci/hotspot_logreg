task merge {
    Array[File]+ vcf_files
    Array[File]+ vcf_files_tbi
    File? use_header
    File? regions_file
    Boolean force_samples = true
    Boolean print_header = false
    Boolean no_version = false
    String? merge
    String output_type = 'v'
    String output_file_prefix
    String outputDir
    Int threads = 3

    command {
        bcftools merge \
        ${true="--force-samples" false="" force_samples} \
        ${true="--print-header" false="" print_header} \
        ${true="--no-version" false="" no_version} \
        ${"--regions-file " + regions_file} \
        ${"--use-header " + use_header} \
        ${"--merge " + merge} \
        ${"--threads " + threads} \
        --output-type ${output_type} \
        --output ${outputDir}${output_file_prefix}.vcf \
        ${sep=" " vcf_files}
    }

    output {
        File merged_vcf = "${output_file_prefix}.vcf"
    }

    runtime {
        tool: "bcftools"
	strategy: "localized"
    }

}

task process_vcf {
     File vcfFile
     String outputFilePrefix
     String outputDir 

     command <<<
         /bin/bash -c 'cat ${vcfFile} | awk '"'"'{print $1"\t"$2"\t"$2"\t"$4"\t"$5}'"'"' > ${outputDir}${outputFilePrefix}.txt'
     >>>

     output {
         File processedVcf = "${outputFilePrefix}.txt"
     }

     runtime {
         tool: "annovar"
         strategy: "localized"
     }
}

task table_annovar {
   File annovarInput
   File annovarDatabase
   Array[String]+ protocol = ["ensGene", "dbnsfp30a"]
   Array[String]+ operation = ["g", "f"]
   String build = "hg19"
   String nastring = "."
   String outputFilePrefix = "table_annovar"
   String outputDir 

   command {
       table_annovar.pl ${annovarInput} \
                        ${annovarDatabase} \
                        --protocol ${sep="," protocol} \
                        --operation ${sep="," operation} \
                        --build ${build} \
                        --nastring ${nastring} \
                        --outfile ${outputDir}${outputFilePrefix}
   }

   output {
       File multianno = "${outputFilePrefix}.${build}_multianno.txt"
   }

   runtime {
	tool: "annovar"
	strategy: "localized"
   }
}

task oncodrivefm_input_gen {
     File vcfFile
     File multiannoFile
     Int vcfHeader = 30
     Int multiannoHeader = 0
     String outputFilePrefix = "oncodrive_inputgen"
     String outputDir 

    command <<<
         /bin/bash -c 'python <<CODE

         import pandas as pd

         # read files
         vcf = pd.read_table("${vcfFile}", header=int(${vcfHeader}), na_values="./.:.:.", low_memory=False)
         multianno = pd.read_table("${multiannoFile}", header=int(${multiannoHeader}), comment="#")

         # drop extra columns
         vcf.drop(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], axis=1, inplace=True)
         multianno = multianno[["Gene.ensGene", "SIFT_score", "Polyphen2_HDIV_score", "MutationAssessor_score"]]

         # merge files
         merged = pd.concat([vcf, multianno], axis=1)

         # reshape
         melted = pd.melt(merged, id_vars=["Gene.ensGene", "SIFT_score", "Polyphen2_HDIV_score", "MutationAssessor_score"], var_name = "sample", value_name = "mutation_status")

         # drop na, reorder columns and rename columns
         filtered = melted.dropna()
         reordered = filtered[["sample", "Gene.ensGene", "SIFT_score", "Polyphen2_HDIV_score", "MutationAssessor_score"]]
         reordered.columns = ["SAMPLE", "GENE", "SIFT", "PPH2", "MA"]

         # write to output file
         reordered.to_csv("${outputDir}${outputFilePrefix}.tdm", sep="\t", header=True, index=False)

         CODE'
     >>>

     output {
         File oncodrivefm_input = "${outputFilePrefix}.tdm"
     }

     runtime {
         tool: "python-numpy-pandas-scikit-learn"
         strategy: "localized"
     }
}

task chromNamesToEnsembl {
    File inputFile
    String ext = "vcf"
    String outputFilePrefix
    String outputDir 

    command <<<
         /bin/bash -c 'sed '"'"'s/^chr//g'"'"' ${inputFile} > ${outputDir}${outputFilePrefix}_ensemblChroms.${ext}'
    >>>

    output {
        File file_with_ensemblchroms = "${outputFilePrefix}_ensemblChroms.${ext}"
    }

    runtime {
        tool: "vcf2maf"
        strategy: "localized"
    }
}

task vcf2maf {
    File inputVCF
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai
    String ncbiBuild
    String outputDir 
    String outputFilePrefix

    command {
        perl /home/vcf2maf.pl \
        --input-vcf ${inputVCF} \
        --output-maf ${outputDir}${outputFilePrefix}.maf \
        --ref-fasta ${refFasta} \
        --ncbi-build ${ncbiBuild} \
        --vep-data ${vepOfflineCacheDir}
    }

    output {
        File maf = "${outputFilePrefix}.maf"
    }

    runtime {
        tool: "vcf2maf"
        strategy: "localized"
    }
}

task maf2maf {
    File inputMAF
    File vepOfflineCacheDir
    File refFasta
    File refFastaFai
    String ncbiBuild
    String outputDir 
    String outputFilePrefix = "minimal_maf"

    command {
        perl /home/maf2maf.pl \
        --input-maf ${inputMAF} \
        --output-maf ${outputDir}${outputFilePrefix}.maf \
        --ref-fasta ${refFasta} \
        --ncbi-build ${ncbiBuild} \
        --vep-data ${vepOfflineCacheDir}
    }

    output {
        File maf = "${outputFilePrefix}.maf"
    }

    runtime {
        tool: "vcf2maf"
        strategy: "localized"
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
     String naStringRegEx = "\\."
     String outputDir 
     String outputFilePrefix = "taylor"

     command <<<
         /bin/bash -c ' \

         # link ref files to expected location \
         ln -s ${refFasta} /mnt/; \
         ln -s ${refFastaFai} /mnt/; \

         # replace NA strings with what the algorithm expects and remove chr \
         sed '"'"'s/${naStringRegEx}/NA/g'"'"' ${inputMAF} | sed '"'"'s#chr##'"'"' > tmp.maf; \

         # call hotspot algorithm \
         Rscript /home/hotspot_algo.R --input-maf=tmp.maf \
                                --rdata=/home/hotspot_algo.Rdata \
                                ${"--gene-query=" + geneQuery} \
                                ${"--align100mer=" + align100mer} \
                                ${"--align24mer=" + align24mer} \
                                --homopolymer=${homopolymer} \
                                --filter-centerbias=${filterCenterBias} \
                                --output-file=${outputDir}${outputFilePrefix}.txt;

        # cleanup \
        rm tmp.maf /mnt/$(basename ${refFasta}) /mnt/$(basename ${refFastaFai})'
     >>>

     output {
         File hotspots = "${outputFilePrefix}.txt"
     }

     runtime {
         tool: "taylor-lab-hotspots"
         strategy: "localized"
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
    String outputDir 

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
        -o ${outputDir} \
        ${inputData}
    }

    output {
        File genes = "${outputFilePrefix}-genes.tsv"
        File pathways = "${outputFilePrefix}-pathways.tsv"
    }

    runtime {
	tool: "oncodrivefm"
	strategy: "localized"
    }
}


workflow hotspot {

    # OutputDir argument from JSON file that \
    # will be passed commonly to all tasks
    String outputDir

    # merge the VCF files
    call merge {
        input: outputDir = outputDir
    }

    ## ONCODRIVE ##

    call process_vcf {
        input: vcfFile = merge.merged_vcf,
               outputDir = outputDir
    }

    call table_annovar {
        input: annovarInput = process_vcf.processedVcf,
               outputDir = outputDir
    }

    call oncodrivefm_input_gen {
        input: vcfFile = merge.merged_vcf,
               multiannoFile = table_annovar.multianno,
               outputDir = outputDir
    }

    call oncodrivefm {
        input: inputData = oncodrivefm_input_gen.oncodrivefm_input, 
               outputDir = outputDir
    }

    ## Taylor ##
    call chromNamesToEnsembl {
        input:  inputFile = merge.merged_vcf,
                outputDir = outputDir
    }

    call vcf2maf {
        input:  inputVCF = chromNamesToEnsembl.file_with_ensemblchroms,
                outputDir = outputDir
    }

    call maf2maf {
        input: inputMAF = vcf2maf.maf,
               outputDir = outputDir
    }

    call taylor {
        input: inputMAF = maf2maf.maf,
               outputDir = outputDir
    }
}
