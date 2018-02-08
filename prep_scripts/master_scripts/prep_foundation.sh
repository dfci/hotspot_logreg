#!/usr/bin/env bash

DATA_DIR='/data/zwiesler/intelccc/data/public/foundation'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- foundation_one_breast_tumors.txt'
echo '- T5a_genes.txt'

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical and genomic data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_foundation.py \
    --in-clinical ${DATA_DIR}/foundation_one_breast_tumors.txt \
    --out-clinical ${DATA_DIR}/foundation.clinical.r1.tsv \
    --out-genomic ${DATA_DIR}/foundation.r1.maf

echo; echo;
echo '-------------------------------'
echo 'Create bed file'
echo '-------------------------------'
${SCRIPT_DIR}/prep_foundation_bed.py \
    --in-bed ${DATA_DIR}/genie_combined.bed \
    --in-genelist ${DATA_DIR}/T5a_genes.txt \
    --out-bed ${DATA_DIR}/foundation.r1.bed

echo; echo;
echo '-------------------------------'
echo 'Validating clinical, genomic, and .bed files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/foundation.clinical.r1.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/foundation.r1.maf
${VALIDATE_DIR}/validate_bed.py --in_bed ${DATA_DIR}/foundataion.r1.bed

echo; echo;
echo '## -------------------------------'
echo '## clinical: ' ${DATA_DIR}/foundation.clinical.r1.tsv
echo '## genomic: ' ${DATA_DIR}/foundation.r1.maf
echo '## -------------------------------'