#!/bin/bash -eu


echo <<EOF

SYNOPSIS

    This bash script perform post-processing on the Mutation Hotspot Table. It
    calls various python subscripts which add new annotation columns to the
    Mutation Hotspot Table.

EXAMPLES

    ./run_postproc.sh \
        ../../private/reliance_point_outputs/05May2017/05May2017_ReliancePoint_ClinicalHotspotsTable.txt \
        ../../private/count_tables/05May2017.count_table.merged.only_public.with_breakdown.no_qvalue_filter.min_count_1.txt \
        ../../private/tmp/final.postproc.05May2017_ReliancePoint_ClinicalHotspotsTable.txt

AUTHOR
    Parin Sripakdeevong <parin@jimmy.harvard.edu> (May-2017)
EOF

if [[ $# -ne 3 ]] ; then

    echo 'USAGE: run_post_proc.sh <in_hotspot_table> <in_public_hotspot_table> <out_hotspot_table>'
    exit 1
fi

echo ">> Enter run_proc.sh."


# Input Options
in_hotspot_table=$1
in_public_hotspot_table=$2
out_hotspot_table=$3

# Inferred Variables
workdir=`dirname $out_hotspot_table`

tmp1_hotspot_table=${workdir}/"tmp1.postproc.`basename ${in_hotspot_table}`"
tmp2_hotspot_table=${workdir}/"tmp2.postproc.`basename ${in_hotspot_table}`"

echo ">> -------------------------------------------"
echo ">> Input_Options:"
echo ">>   in_hotspot_table ${in_hotspot_table}"
echo ">>   in_public_hotspot_table: ${in_public_hotspot_table}"
echo ">>   out_hotspot_table: ${out_hotspot_table}"
echo ">> -------------------------------------------"
echo ">> Derived_Variables:"
echo ">>   workdir ${workdir}"
echo ">>   tmp1_hotspot_table: ${tmp1_hotspot_table}"
echo ">>   tmp2_hotspot_table: ${tmp2_hotspot_table}"
echo ">> -------------------------------------------"

echo ""
echo ">> Create workdir"
mkdir -v -p ${workdir}

echo ""
echo ">> Run add_taylor_cols.py"
./add_taylor_cols.py \
    --in_table ${in_hotspot_table} \
    --out_table ${tmp1_hotspot_table}

echo ""
echo ">> Run add_clin_qval_cols.py"
./add_clin_qval_cols.py \
        --in_table ${tmp1_hotspot_table} \
        --out_table ${tmp2_hotspot_table}

echo ""
echo ">> Run add_public_only_cols.py"
 ./add_public_only_cols.py \
        --public_only_table ${in_public_hotspot_table} \
        --in_table ${tmp2_hotspot_table} \
        --out_table ${out_hotspot_table}

echo ""
echo ">> Exit run_postproc.sh."
