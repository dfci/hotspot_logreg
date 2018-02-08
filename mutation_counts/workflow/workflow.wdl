workflow wf {
  Array[File]+ inMafs
  Array[File]+ inTargetBeds
  Array[File]+ inClinicals
  File scriptsDir
  File geneAnnotationGtf
  File refFasta

  call getSourceMafIDs {
    input:
      inMafs=inMafs
  }

  call mergeAnnotatedMaf {
    input:
      inMafs=inMafs,
      scriptsDir=scriptsDir
  }

  call mergeClinical {
    input:
      inClinicals=inClinicals,
      scriptsDir=scriptsDir
  }

  call countMutations {
    input:
      inMaf=mergeAnnotatedMaf.out,
      inClinical=mergeClinical.out,
      inSourceMafIDs=getSourceMafIDs.out,
      inTargetBeds=inTargetBeds,
      scriptsDir=scriptsDir,
      geneAnnotationGtf=geneAnnotationGtf,
      refFasta=refFasta,
      refFastaFai=refFasta + ".fai"
  }
}

task mergeAnnotatedMaf {
  Array[File]+ inMafs
  File scriptsDir

  command {
    ${scriptsDir}/merge_annotated_maf.py \
      --in_mafs ${sep=" " inMafs} \
      --out_maf out
  }

  output {
    File out = "out"
  }

  runtime {
    docker: "dfciksg/pipeline-python:v1.0.6"
  }
}

task mergeClinical {
  Array[File]+ inClinicals
  File scriptsDir

  command {
    ${scriptsDir}/merge_clinical.py \
      --in_clinicals ${sep=" " inClinicals} \
      --out_clinical out
  }

  output {
    File out = "out"
  }

  runtime {
    docker: "dfciksg/pipeline-python:v1.0.6"
  }
}

task getSourceMafIDs {
  Array[File]+ inMafs

  command <<<
    python <<CODE
    import os
    in_maf_files = "${sep="," inMafs}".split(",")
    for in_maf_file in in_maf_files:
        print os.path.basename(in_maf_file)
    CODE
  >>>

  output {
    Array[String] out = read_lines(stdout())
  }

  runtime {
    docker: "dfciksg/pipeline-python:v1.0.6"
  }
}

task countMutations {
  File inMaf
  File inClinical
  Array[String]+ inSourceMafIDs
  Array[File]+ inTargetBeds
  File scriptsDir
  File geneAnnotationGtf
  File refFasta
  File refFastaFai
  Int minCount
  Boolean runHotspotStatsModel
  Boolean useToyStatsModel
  String biopsySiteTypeFilter

  command {
    ${scriptsDir}/count_mutations.py \
      --in_maf ${inMaf} \
      --in_clinical ${inClinical} \
      --in_source_maf_ids ${sep=" " inSourceMafIDs} \
      --in_target_beds ${sep=" " inTargetBeds} \
      --gene_annotation_gtf ${geneAnnotationGtf} \
      --ref_fasta ${refFasta} \
      --ref_fasta_fai ${refFastaFai} \
      --min_count ${minCount} \
      ${true="--run_hotspot_stats_model" false="" runHotspotStatsModel} \
      ${true="--use_toy_stats_model" false="" useToyStatsModel} \
      --biopsy_site_type_filter ${biopsySiteTypeFilter} \
      --out_table count_table.txt
  }

  output {
    File out = "count_table.txt"
  }

  runtime {
    docker: "dfciksg/pipeline-python:v1.0.6"
  }
}
