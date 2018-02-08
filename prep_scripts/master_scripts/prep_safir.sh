#!/usr/bin/env bash

DATA_DIR='/data/zwiesler/intelccc/data/public/safir'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- data_clinical_patient.txt'
echo '- data_clinical_sample.txt'
echo '- data_mutations_extended.txt'

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_safir_clinical.py \
    --in-sample ${DATA_DIR}/data_clinical_sample.txt \
    --out-clinical ${DATA_DIR}/safir.clinical.r1.tsv

echo; echo;
echo '-------------------------------'
echo 'Normalize genomic data'
echo '-------------------------------'
echo '## NOTE: I am using prep_tcga_maf even though this '
echo '##       script was build for a different dataset '
echo '##       because both datasets were downloaded from '
echo '##       cBioPortal and the .maf formats should be consistent.'; echo;
${SCRIPT_DIR}/prep_safir_genomic.py \
    --in-maf ${DATA_DIR}/data_mutations_extended.txt \
    --out-genomic ${DATA_DIR}/safir.r1.maf

echo; echo;
echo '-------------------------------'
echo 'Validating clinical and genomic files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/safir.clinical.r1.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/safir.r1.maf

echo; echo;
echo '## -------------------------------'
echo '## clinical: ' ${DATA_DIR}/safir.clinical.r1.tsv
echo '## genomic: ' ${DATA_DIR}/safir.r1.maf
echo '## -------------------------------'