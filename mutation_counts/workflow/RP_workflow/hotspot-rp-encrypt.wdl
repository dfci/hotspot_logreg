task DFCI_encrypt_file {
  File DFCI_input_maf_file_1
  File DFCI_input_maf_file_2
  File DFCI_input_maf_file_3
  File DFCI_input_bed_file_1
  File DFCI_input_bed_file_2
  File DFCI_input_bed_file_3
  File DFCI_input_tsv_file_1
  String DFCI_encrypted_maf_output_1
  String DFCI_encrypted_maf_output_2
  String DFCI_encrypted_maf_output_3
  String DFCI_encrypted_bed_output_1
  String DFCI_encrypted_bed_output_2
  String DFCI_encrypted_bed_output_3
  String DFCI_encrypted_tsv_output_1

  command {
    ./encrypt.py -i ${DFCI_input_maf_file_1} ${DFCI_input_maf_file_2} ${DFCI_input_maf_file_3} ${DFCI_input_bed_file_1} ${DFCI_input_bed_file_2} ${DFCI_input_bed_file_3} ${DFCI_input_tsv_file_1} -o /output/${DFCI_encrypted_maf_output_1} /output/${DFCI_encrypted_maf_output_2} /output/${DFCI_encrypted_maf_output_3} /output/${DFCI_encrypted_bed_output_1} /output/${DFCI_encrypted_bed_output_2} /output/${DFCI_encrypted_bed_output_3} /output/${DFCI_encrypted_tsv_output_1}
  }

  output {
    File DFCI_response_maf_1 = "${DFCI_encrypted_maf_output_1}"
    File DFCI_response_maf_2 = "${DFCI_encrypted_maf_output_2}"
    File DFCI_response_maf_3 = "${DFCI_encrypted_maf_output_3}"
    File DFCI_response_bed_1 = "${DFCI_encrypted_bed_output_1}"
    File DFCI_response_bed_2 = "${DFCI_encrypted_bed_output_2}"
    File DFCI_response_bed_3 = "${DFCI_encrypted_bed_output_3}"
    File DFCI_response_tsv_1 = "${DFCI_encrypted_tsv_output_1}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task DFCI_encrypt_public_file1 {
  File DFCI_IMPACT341_maf
  File DFCI_IMPACT410_maf
  File DFCI_IMPACT341_bed
  File DFCI_IMPACT410_bed
  File DFCI_IMPACT341_tsv
  File DFCI_IMPACT410_tsv
  String DFCI_IMPACT341_encrypted_maf
  String DFCI_IMPACT410_encrypted_maf
  String DFCI_IMPACT341_encrypted_bed
  String DFCI_IMPACT410_encrypted_bed
  String DFCI_IMPACT341_encrypted_tsv
  String DFCI_IMPACT410_encrypted_tsv
  
  command {
    ./encrypt.py -i ${DFCI_IMPACT341_maf} ${DFCI_IMPACT410_maf} ${DFCI_IMPACT341_bed} ${DFCI_IMPACT410_bed} ${DFCI_IMPACT341_tsv} ${DFCI_IMPACT410_tsv} -o /output/${DFCI_IMPACT341_encrypted_maf} /output/${DFCI_IMPACT410_encrypted_maf} /output/${DFCI_IMPACT341_encrypted_bed} /output/${DFCI_IMPACT410_encrypted_bed} /output/${DFCI_IMPACT341_encrypted_tsv} /output/${DFCI_IMPACT410_encrypted_tsv}
  }

  output {                                                      
    File DFCI_response_maf_1 = "${DFCI_IMPACT341_encrypted_maf}" 
    File DFCI_response_maf_2 = "${DFCI_IMPACT410_encrypted_maf}" 
    File DFCI_response_bed_1 = "${DFCI_IMPACT341_encrypted_bed}" 
    File DFCI_response_bed_2 = "${DFCI_IMPACT410_encrypted_bed}" 
    File DFCI_response_tsv_1 = "${DFCI_IMPACT341_encrypted_tsv}"
	File DFCI_response_tsv_2 = "${DFCI_IMPACT410_encrypted_tsv}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task DFCI_encrypt_public_file2 {
  File DFCI_METABRIC_maf
  File DFCI_METABRIC_bed
  File DFCI_METABRIC_tsv
  String DFCI_METABRIC_encrypted_maf
  String DFCI_METABRIC_encrypted_bed
  String DFCI_METABRIC_encrypted_tsv
  File DFCI_TCGA_BRCA_maf
  File DFCI_TCGA_BRCA_bed
  File DFCI_TCGA_BRCA_tsv
  String DFCI_TCGA_BRCA_encrypted_maf
  String DFCI_TCGA_BRCA_encrypted_bed
  String DFCI_TCGA_BRCA_encrypted_tsv
  File DFCI_Sanger_maf
  File DFCI_Sanger_bed
  File DFCI_Sanger_tsv
  String DFCI_Sanger_encrypted_maf
  String DFCI_Sanger_encrypted_bed
  String DFCI_Sanger_encrypted_tsv


  command {
    ./encrypt.py -i ${DFCI_METABRIC_maf} ${DFCI_METABRIC_bed} ${DFCI_METABRIC_tsv} ${DFCI_TCGA_BRCA_maf} ${DFCI_TCGA_BRCA_bed} ${DFCI_TCGA_BRCA_tsv} ${DFCI_Sanger_maf} ${DFCI_Sanger_bed} ${DFCI_Sanger_tsv} -o /output/${DFCI_METABRIC_encrypted_maf} /output/${DFCI_METABRIC_encrypted_bed} /output/${DFCI_METABRIC_encrypted_tsv} /output/${DFCI_TCGA_BRCA_encrypted_maf} /output/${DFCI_TCGA_BRCA_encrypted_bed} /output/${DFCI_TCGA_BRCA_encrypted_tsv} /output/${DFCI_Sanger_encrypted_maf} /output/${DFCI_Sanger_encrypted_bed} /output/${DFCI_Sanger_encrypted_tsv}
  }

  output {
    File DFCI_METABRIC_maf_response = "${DFCI_METABRIC_encrypted_maf}"
    File DFCI_METABRIC_bed_response = "${DFCI_METABRIC_encrypted_bed}"
    File DFCI_METABRIC_tsv_response = "${DFCI_METABRIC_encrypted_tsv}"
    File DFCI_TCGA_BRCA_maf_response = "${DFCI_TCGA_BRCA_encrypted_maf}"
    File DFCI_TCGA_BRCA_bed_response = "${DFCI_TCGA_BRCA_encrypted_bed}"
    File DFCI_TCGA_BRCA_tsv_response = "${DFCI_TCGA_BRCA_encrypted_tsv}"
    File DFCI_Sanger_maf_response = "${DFCI_Sanger_encrypted_maf}"
    File DFCI_Sanger_bed_response = "${DFCI_Sanger_encrypted_bed}"
    File DFCI_Sanger_tsv_response = "${DFCI_Sanger_encrypted_tsv}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task DFCI_encrypt_public_file3 {
  File DFCI_icgc_brca_uk_vep_maf
  File DFCI_icgc_brca_fr_vep_maf
  File DFCI_icgc_brca_kr_vep_maf
  File DFCI_safir_r1_vep_maf
  File DFCI_genie_grcc_r1_vep_maf
  File DFCI_genie_vicc_r1_vep_maf
  File DFCI_genie_mda_r1_vep_maf
  File DFCI_genie_uhn_r1_vep_maf
  File DFCI_foundation_r1_vep_maf
  File DFCI_GENIE_grcc_r1_bed
  File DFCI_GENIE_vicc_r1_bed
  File DFCI_GENIE_mda_r1_bed
  File DFCI_GENIE_uhn_r1_bed
  File DFCI_foundation_r1_bed
  File DFCI_icgc_brca_uk_clinical_r1_tsv
  File DFCI_icgc_brca_fr_clinical_r1_tsv
  File DFCI_icgc_brca_kr_clinical_r1_tsv
  File DFCI_safir_clinical_r1_tsv
  File DFCI_genie_grcc_clinical_r1_tsv
  File DFCI_genie_vicc_clinical_r1_tsv
  File DFCI_genie_mda_clinical_r1_tsv
  File DFCI_genie_uhn_clinical_r1_tsv
  File DFCI_foundation_clinical_r1_tsv
  String DFCI_encrypted_icgc_brca_uk_vep_maf
  String DFCI_encrypted_icgc_brca_fr_vep_maf
  String DFCI_encrypted_icgc_brca_kr_vep_maf
  String DFCI_encrypted_safir_r1_vep_maf
  String DFCI_encrypted_genie_grcc_r1_vep_maf
  String DFCI_encrypted_genie_vicc_r1_vep_maf
  String DFCI_encrypted_genie_mda_r1_vep_maf
  String DFCI_encrypted_genie_uhn_r1_vep_maf
  String DFCI_encrypted_foundation_r1_vep_maf
  String DFCI_encrypted_GENIE_grcc_r1_bed
  String DFCI_encrypted_GENIE_vicc_r1_bed
  String DFCI_encrypted_GENIE_mda_r1_bed
  String DFCI_encrypted_GENIE_uhn_r1_bed
  String DFCI_encrypted_foundation_r1_bed
  String DFCI_encrypted_icgc_brca_uk_clinical_r1_tsv
  String DFCI_encrypted_icgc_brca_fr_clinical_r1_tsv
  String DFCI_encrypted_icgc_brca_kr_clinical_r1_tsv
  String DFCI_encrypted_safir_clinical_r1_tsv
  String DFCI_encrypted_genie_grcc_clinical_r1_tsv
  String DFCI_encrypted_genie_vicc_clinical_r1_tsv
  String DFCI_encrypted_genie_mda_clinical_r1_tsv
  String DFCI_encrypted_genie_uhn_clinical_r1_tsv
  String DFCI_encrypted_foundation_clinical_r1_tsv

  command {
    ./encrypt.py -i ${DFCI_icgc_brca_uk_vep_maf} ${DFCI_icgc_brca_fr_vep_maf} ${DFCI_icgc_brca_kr_vep_maf} ${DFCI_safir_r1_vep_maf} ${DFCI_genie_grcc_r1_vep_maf} ${DFCI_genie_vicc_r1_vep_maf} ${DFCI_genie_mda_r1_vep_maf} ${DFCI_genie_uhn_r1_vep_maf} ${DFCI_foundation_r1_vep_maf} ${DFCI_GENIE_grcc_r1_bed} ${DFCI_GENIE_vicc_r1_bed} ${DFCI_GENIE_mda_r1_bed} ${DFCI_GENIE_uhn_r1_bed} ${DFCI_foundation_r1_bed} ${DFCI_icgc_brca_uk_clinical_r1_tsv} ${DFCI_icgc_brca_fr_clinical_r1_tsv} ${DFCI_icgc_brca_kr_clinical_r1_tsv} ${DFCI_safir_clinical_r1_tsv} ${DFCI_genie_grcc_clinical_r1_tsv} ${DFCI_genie_vicc_clinical_r1_tsv} ${DFCI_genie_mda_clinical_r1_tsv} ${DFCI_genie_uhn_clinical_r1_tsv} ${DFCI_foundation_clinical_r1_tsv} -o /output/${DFCI_encrypted_icgc_brca_uk_vep_maf} /output/${DFCI_encrypted_icgc_brca_fr_vep_maf} /output/${DFCI_encrypted_icgc_brca_kr_vep_maf} /output/${DFCI_encrypted_safir_r1_vep_maf} /output/${DFCI_encrypted_genie_grcc_r1_vep_maf} /output/${DFCI_encrypted_genie_vicc_r1_vep_maf} /output/${DFCI_encrypted_genie_mda_r1_vep_maf} /output/${DFCI_encrypted_genie_uhn_r1_vep_maf} /output/${DFCI_encrypted_foundation_r1_vep_maf} /output/${DFCI_encrypted_GENIE_grcc_r1_bed} /output/${DFCI_encrypted_GENIE_vicc_r1_bed} /output/${DFCI_encrypted_GENIE_mda_r1_bed} /output/${DFCI_encrypted_GENIE_uhn_r1_bed} /output/${DFCI_encrypted_foundation_r1_bed} /output/${DFCI_encrypted_icgc_brca_uk_clinical_r1_tsv} /output/${DFCI_encrypted_icgc_brca_fr_clinical_r1_tsv} /output/${DFCI_encrypted_icgc_brca_kr_clinical_r1_tsv} /output/${DFCI_encrypted_safir_clinical_r1_tsv} /output/${DFCI_encrypted_genie_grcc_clinical_r1_tsv} /output/${DFCI_encrypted_genie_vicc_clinical_r1_tsv} /output/${DFCI_encrypted_genie_mda_clinical_r1_tsv} /output/${DFCI_encrypted_genie_uhn_clinical_r1_tsv} /output/${DFCI_encrypted_foundation_clinical_r1_tsv}
  }

  output {
   File response1 = "${DFCI_encrypted_icgc_brca_uk_vep_maf}"
   File response2 = "${DFCI_encrypted_icgc_brca_fr_vep_maf}"
   File response3 = "${DFCI_encrypted_icgc_brca_kr_vep_maf}"
   File response4 = "${DFCI_encrypted_safir_r1_vep_maf}"
   File response5 = "${DFCI_encrypted_genie_grcc_r1_vep_maf}"
   File response6 = "${DFCI_encrypted_genie_vicc_r1_vep_maf}"
   File response7 = "${DFCI_encrypted_genie_mda_r1_vep_maf}"
   File response8 = "${DFCI_encrypted_genie_uhn_r1_vep_maf}"
   File response9 = "${DFCI_encrypted_foundation_r1_vep_maf}"
   File response10 = "${DFCI_encrypted_GENIE_grcc_r1_bed}"
   File response11 = "${DFCI_encrypted_GENIE_vicc_r1_bed}"
   File response12 = "${DFCI_encrypted_GENIE_mda_r1_bed}"
   File response13 = "${DFCI_encrypted_GENIE_uhn_r1_bed}"
   File response14 = "${DFCI_encrypted_foundation_r1_bed}"
   File response15 = "${DFCI_encrypted_icgc_brca_uk_clinical_r1_tsv}"
   File response16 = "${DFCI_encrypted_icgc_brca_fr_clinical_r1_tsv}"
   File response17 = "${DFCI_encrypted_icgc_brca_kr_clinical_r1_tsv}"
   File response18 = "${DFCI_encrypted_safir_clinical_r1_tsv}"
   File response19 = "${DFCI_encrypted_genie_grcc_clinical_r1_tsv}"
   File response20 = "${DFCI_encrypted_genie_vicc_clinical_r1_tsv}"
   File response21 = "${DFCI_encrypted_genie_mda_clinical_r1_tsv}"
   File response22 = "${DFCI_encrypted_genie_uhn_clinical_r1_tsv}"
   File response23 = "${DFCI_encrypted_foundation_clinical_r1_tsv}"
 }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task DFCI_encrypt_resource_file {
  File DFCI_fasta
  File DFCI_fai
  File DFCI_gtf
  String DFCI_encrypted_fasta
  String DFCI_encrypted_fai
  String DFCI_encrypted_gtf

  command {
    ./encrypt.py -i ${DFCI_fasta} ${DFCI_fai} ${DFCI_gtf} -o /output/${DFCI_encrypted_fasta} /output/${DFCI_encrypted_fai} /output/${DFCI_encrypted_gtf}
  }

  output {
    File DFCI_response_fasta = "${DFCI_encrypted_fasta}"
    File DFCI_response_fai = "${DFCI_encrypted_fai}"
    File DFCI_response_gtf = "${DFCI_encrypted_gtf}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task OICR_encrypt_file {
  File OICR_input_maf_file_1
  File OICR_input_bed_file_1
  File OICR_input_tsv_file_1
  String OICR_encrypted_maf_output_1
  String OICR_encrypted_bed_output_1
  String OICR_encrypted_tsv_output_1

  command {
    ./encrypt.py -i ${OICR_input_maf_file_1} ${OICR_input_bed_file_1} ${OICR_input_tsv_file_1} -o /output/${OICR_encrypted_maf_output_1} /output/${OICR_encrypted_bed_output_1} /output/${OICR_encrypted_tsv_output_1}
  }

  output {
    File OICR_response_maf_1 = "${OICR_encrypted_maf_output_1}"
    File OICR_response_bed_1 = "${OICR_encrypted_bed_output_1}"
    File OICR_response_tsv_1 = "${OICR_encrypted_tsv_output_1}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

task OHSU_encrypt_file {
  File OHSU_input_maf_file_1
  File OHSU_input_bed_file_1
  File OHSU_input_tsv_file_1
  String OHSU_encrypted_maf_output_1
  String OHSU_encrypted_bed_output_1
  String OHSU_encrypted_tsv_output_1


  command {
    ./encrypt.py -i ${OHSU_input_maf_file_1} ${OHSU_input_bed_file_1} ${OHSU_input_tsv_file_1}  -o /output/${OHSU_encrypted_maf_output_1} /output/${OHSU_encrypted_bed_output_1} /output/${OHSU_encrypted_tsv_output_1} 
  }

  output {
    File OHSU_response_maf_1 = "${OHSU_encrypted_maf_output_1}"
    File OHSU_response_bed_1 = "${OHSU_encrypted_bed_output_1}"
    File OHSU_response_tsv_1 = "${OHSU_encrypted_tsv_output_1}"
  }

  runtime {
    tool: "rp_encrypt_app"
    strategy: "fetch_file"
  }
}

workflow encrypt_hotspots {
  call DFCI_encrypt_file
  call DFCI_encrypt_public_file1
  call DFCI_encrypt_public_file2
  call DFCI_encrypt_public_file3
  call DFCI_encrypt_resource_file
  call OICR_encrypt_file
  call OHSU_encrypt_file
}
