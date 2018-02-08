#!/usr/bin/env bash

DATA_DIR='/data/zwiesler/intelccc/data/public/genie/genie_other'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- data_clinical.txt'
echo '- data_mutations_extended.txt'

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_genie_clinical.py \
    --in-clinical ${DATA_DIR}/data_clinical.txt \
    --out-clinical ${DATA_DIR}

echo; echo;
echo '-------------------------------'
echo 'Normalize genomic data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_genie_genomic.py \
    --in-maf ${DATA_DIR}/data_mutations_extended.txt \
    --in-clinical ${DATA_DIR} \
    --out-genomic ${DATA_DIR}

echo; echo;
echo '-------------------------------'
echo 'Create bed file'
echo '-------------------------------'
${SCRIPT_DIR}/prep_genie_bed.py \
    --in-bed ${DATA_DIR}/genie_combined.bed \
    --out-dir ${DATA_DIR}

echo; echo;
echo '-------------------------------'
echo 'Validating clinical, genomic, and bed files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/GRCC.clinical.r1.tsv
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/VICC.clinical.r1.tsv
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/MDA.clinical.r1.tsv
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/UHN.clinical.r1.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/GENIE.grcc.r1.maf
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/GENIE.vicc.r1.maf
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/GENIE.mda.r1.maf
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/GENIE.uhn.r1.maf
${VALIDATE_DIR}/validate_bed.py --in_bed ${DATA_DIR}/GENIE.grcc.r1.bed
${VALIDATE_DIR}/validate_bed.py --in_bed ${DATA_DIR}/GENIE.vicc.r1.bed
${VALIDATE_DIR}/validate_bed.py --in_bed ${DATA_DIR}/GENIE.mda.r1.bed
${VALIDATE_DIR}/validate_bed.py --in_bed ${DATA_DIR}/GENIE.uhn.r1.bed

echo; echo;
echo '## -------------------------------'
echo '## clinical: ' ${DATA_DIR}/GRCC.clinical.r1.tsv
echo '## clinical: ' ${DATA_DIR}/VICC.clinical.r1.tsv
echo '## clinical: ' ${DATA_DIR}/MDA.clinical.r1.tsv
echo '## clinical: ' ${DATA_DIR}/UHN.clinical.r1.tsv
echo '## genomic: ' ${DATA_DIR}/GENIE.grcc.r1.maf
echo '## genomic: ' ${DATA_DIR}/GENIE.vicc.r1.maf
echo '## genomic: ' ${DATA_DIR}/GENIE.mda.r1.maf
echo '## genomic: ' ${DATA_DIR}/GENIE.uhn.r1.maf
echo '## -------------------------------'