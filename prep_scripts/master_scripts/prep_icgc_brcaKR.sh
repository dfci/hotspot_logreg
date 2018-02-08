#!/usr/bin/env bash

DATA_DIR='/data/zwiesler/intelccc/data/public/icgc/brca-kr'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- donor.BRCA-KR.tsv'
echo '- sample.BRCA-KR.tsv'
echo '- simple_somatic_mutation.open.BRCA-KR.tsv'

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_icgc_brcaKR_clinical.py \
    --donor-file ${DATA_DIR}/donor.BRCA-KR.tsv \
    --sample-file ${DATA_DIR}/sample.BRCA-KR.tsv \
    --out-clinical ${DATA_DIR}/icgc_brca_kr.clinical.r1.tsv

echo; echo;
echo '-------------------------------'
echo 'Normalize genomic data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_icgc_genomic.py \
    --in-vcf ${DATA_DIR}/simple_somatic_mutation.open.BRCA-KR.tsv \
    --out-maf ${DATA_DIR}/icgc_brca_kr.r1.vep.maf

echo; echo;
echo '-------------------------------'
echo 'Validating clinical and genomic files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/icgc_brca_kr.clinical.r1.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/icgc_brca_kr.r1.vep.maf

echo; echo;
echo '## -------------------------------'
echo '## clinical: '${DATA_DIR}/icgc_brca_kr.clinical.r1.tsv
echo '## genomic: ' ${DATA_DIR}/icgc_brca_kr.r1.vep.maf
echo '## -------------------------------'