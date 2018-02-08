#!/usr/bin/env bash

INSIGHT_DIR='/data/pipeline/v1/modules/annotation/latest/bin'
DATA_DIR='/data/zwiesler/intelccc/data/private/22May2017'
SCRIPT_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/prep_scripts'
VALIDATE_DIR='/data/zwiesler/intelccc/CCC_HotspotUsecase/mutation_counts/scripts'

#echo '-------------------------------'
#echo 'Pulling CAMD data from CRDR API'
#echo '-------------------------------'
## TODO write python script to pull down data
echo 'Data directory is:    ' ${DATA_DIR}
echo 'Script directory is:  ' ${SCRIPT_DIR}
echo; echo 'You must have the following files in your data directory:'
echo '- OP_VARIANT.TXT'
echo '- data_mutations_extended.txt'
echo '- data_samples.txt'
echo '- data_patients.txt'


echo; echo;
echo '-------------------------------'
echo 'Annotate PROFILE data'
echo '-------------------------------'
${INSIGHT_DIR}/annot_profile.py \
    --in_profile ${DATA_DIR}/OP_VARIANT.TXT \
    --in_cbioone ${DATA_DIR}/data_mutations_extended.txt \
    --outdir ./profile_annotated

echo; echo;
echo '-------------------------------'
echo 'Normalize clinical data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_dfci_clinical.py \
    --in_camd_clinical ${DATA_DIR}/data_samples.txt \
    --in_camd_gender ${DATA_DIR}/data_patients.txt \
    --in_haoguo_clinical /data/zwiesler/intelccc/data/private/clinical_de_identified.sample.txt \
    --out_clinical ${DATA_DIR}/dfci.clinical_data.r6.tsv

${SCRIPT_DIR}/prep_dfci_clinical.py \
    --in_camd_clinical ${DATA_DIR}/data_samples.txt \
    --in_camd_gender ${DATA_DIR}/data_patients.txt \
    --in_haoguo_clinical /data/zwiesler/intelccc/data/private/clinical_de_identified.sample.txt \
    --details_mode \
    --out_clinical ${DATA_DIR}/dfci.clinical_data.r6.details.tsv

echo; echo;
echo '-------------------------------'
echo 'Normalize genomic data'
echo '-------------------------------'
${SCRIPT_DIR}/prep_dfci_maf.py \
    --in_genomics ./profile_annotated/vep_annotated.maf \
    --in_clinical ${DATA_DIR}/dfci.clinical_data.r6.details.tsv \
    --panel_versions 1 2 3 \
    --out_mafs ${DATA_DIR}/dfci.DFCI-ONCOPANEL-1.r6.maf \
               ${DATA_DIR}/dfci.DFCI-ONCOPANEL-2.r6.maf \
               ${DATA_DIR}/dfci.DFCI-ONCOPANEL-3.r6.maf

echo; echo;
echo '-------------------------------'
echo 'Validating clinical and genomic files'
echo '-------------------------------'
${VALIDATE_DIR}/validate_clinical.py --in_clinical ${DATA_DIR}/dfci.clinical_data.r5.tsv
${VALIDATE_DIR}/validate_maf.py --in_maf ${DATA_DIR}/dfci.DFCI-ONCOPANEL-1.r6.maf

echo; echo;
echo '## -------------------------------'
echo '## clinical: ' ${DATA_DIR}/dfci.clinical_data.r5.tsv
echo '## genomic: ' ${DATA_DIR}/dfci.DFCI-ONCOPANEL-1.r6..maf
echo '## -------------------------------'
