task mergeAnnotatedMaf {
    Array[File] inMafs
    String outMaf

    command {
        python /opt/rp/submit_job.py -j 1 -r '"-i ${sep="," inMafs} -o ${outMaf}"'
    }

    output {
        File response = outMaf
    }

    runtime {
        tool: "mergeAnnotatedMaf"
	strategy: "enclave"
        arguments: "-i 19"
    }

}

task getSourceMafIDs {
    File Dummy
    Array[File] inMafs
    String sourceMaf_outFile

    command {
        python /opt/rp/submit_job.py -j 2 -r '"-i ${sep="," inMafs} -o ${sourceMaf_outFile}"'
    }
    
     output {
        File response = sourceMaf_outFile
    }

    runtime {
        tool: "getSourceMafIDs"
	strategy: "enclave"
        arguments: "-i 19"
  }
}

task mergeClinical {
    File Dummy
    Array[File] inClinicals
    String outClinical

    command {
        python /opt/rp/submit_job.py -j 3 -r '"-i ${sep="," inClinicals} -o ${outClinical}"'
    }

    output {
        File response = outClinical
    }

    runtime {
        tool: "mergeClinical"
	strategy: "enclave"
        arguments: "-i 17"
    }

}

task countMutations {
    File inMaf
    File inSourceMafIDs
    File inClinical
    Array[File] inTargetBeds
    String mutations_count_output
    File geneAnnotationGtf
    File refFasta
    File refFastaFai

    command {
        python /opt/rp/submit_job.py -j 4 -r '"-i ${inMaf},${inSourceMafIDs},${inClinical},${sep="," inTargetBeds},${geneAnnotationGtf},${refFasta},${refFastaFai} -o ${mutations_count_output}"'
    }

     output {
        File response = mutations_count_output
    }

    runtime {
        tool: "count_mutations"
        strategy: "enclave"
        arguments: "-i 25"
  }
}

workflow RP_Hotspot {
  Array[File] inMafs
  Array[File] inTargetBeds
  Array[File] inClinicals
  File geneAnnotationGtf
  File refFasta
  File refFastaFai

  call mergeAnnotatedMaf { 
    input: 
      	inMafs = inMafs 
  }

  call getSourceMafIDs {
    input: 
       Dummy = mergeAnnotatedMaf.response,
       inMafs = inMafs
  }
  
  call mergeClinical {        
    input:                    
      Dummy = getSourceMafIDs.response,
      inClinicals=inClinicals
  }                           
 
  call countMutations {                       
    input:                                    
      inMaf = mergeAnnotatedMaf.response,       
      inSourceMafIDs = getSourceMafIDs.response,
      inClinical=mergeClinical.response,
      inTargetBeds = inTargetBeds,
      geneAnnotationGtf = geneAnnotationGtf,
      refFasta = refFasta,
      refFastaFai = refFastaFai
  }                                             
}
