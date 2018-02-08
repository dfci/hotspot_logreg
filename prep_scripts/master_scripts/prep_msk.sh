#!/usr/bin/env bash

DATA_DIR='/data/zwiesler/intelccc/data/public/genie/msk_impact_2017'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- data_clinical_patient.txt'
echo '- data_clinical_sample.txt'
echo '- data_mutations_annotated.txt'

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_msk_clinical.py \
    --in-sample ${DATA_DIR}/data_clinical_sample.txt \
    --in-patient ${DATA_DIR}/data_clinical_patient.txt \
    --out-dir ${DATA_DIR}

echo; echo;
echo '-------------------------------'
echo 'Normalize genomic data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_msk_genomic.py \
    --in-maf ${DATA_DIR}/data_mutations_annotated.txt \
    --in-clinical ${DATA_DIR} \
    --out-dir ${DATA_DIR}

echo; echo;
echo '-------------------------------'
echo 'Validating clinical and genomic files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/msk10k.impact341.clinical.r2.tsv
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/msk10k.impact410.clinical.r2.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/msk10k.impact341.r4.maf
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/msk10k.impact410.r4.maf

echo; echo;
echo '## -------------------------------'
echo '## clinical: ' ${DATA_DIR}/msk10k.impact341.clinical.r2.tsv
echo '##         : ' ${DATA_DIR}/msk10k.impact410.clinical.r2.tsv
echo '## genomic: ' ${DATA_DIR}/msk10k.impact341.r4.maf
echo '##        : ' ${DATA_DIR}/msk10k.impact410.r4.maf
echo '## -------------------------------'