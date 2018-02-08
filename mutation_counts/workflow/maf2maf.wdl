workflow wf {
  File inMaf
  File refFasta

  call maf2maf {
    input:
      inFile=inMaf,
      refFasta=refFasta,
      refFastaFai=refFasta + ".fai",
      refFastaGzi=refFasta + ".gzi"
  }
}

task maf2maf {
  File inFile
  File refFasta
  File refFastaFai
  File refFastaGzi
  File vepOfflineCacheTarGz
  String ncbiBuild
  String? customEnst
  String? retainCols

  command {
      mkdir /resource && \
      tar -zxvf ${vepOfflineCacheTarGz} -C /resource && \
      perl /home/maf2maf.pl \
        --input-maf ${inFile} \
        --ref-fasta ${refFasta} \
        --ncbi-build ${ncbiBuild} \
        --vep-data /resource/vep \
        ${'--custom-enst ' + customEnst} \
        ${'--retain-cols ' + retainCols} \
        --output-maf "out.vep.maf"
  }

  output {
    File out = "out.vep.maf"
  }

  runtime {
    docker: "dfciksg/vcf2maf:v1.6.9_VEP85"
  }
}
