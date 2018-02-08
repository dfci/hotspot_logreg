#!/bin/bash -eu

echo ">> Enter pull.sh"
echo ""

google_url="https://storage.googleapis.com"

files=( "vep/v1.0.1_85/vep.tar.gz"
        "ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
        "ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai"
        "ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
        "ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz.fai"
        "ref_genome/GRCh37_r2/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz.gzi"
        "gene_annotation/r1/CDS_only.Homo_sapiens.GRCh37.85.r1.gtf"
        "gene_annotation/r1/subset_genes.CDS_only.Homo_sapiens.GRCh37.85.r1.gtf" )

start_dir=`pwd`

for file_path in "${files[@]}"
do
    dest_dir=`dirname ${file_path}`
    mkdir -v -p ${dest_dir}

    if [ ! -f ${file_path} ]; then
        source="${google_url}/pipeline_resources/${file_path}"
        wget -P ${dest_dir} ${source}
    else
        echo ""
        echo ">> Skip wget. File '${file_path}' already exist."
    fi

    md5_file_path="${file_path}.md5"

    echo ""
    echo ">> Download MD5 file (remove existing)."
    if [ -f ${md5_file_path} ]; then
        rm -r -v ${md5_file_path}
    fi

    source_md5="${google_url}/pipeline_resources/${md5_file_path}"
    wget -P ${dest_dir} ${source_md5}

    echo ""
    echo ">> Perform MD5 check_sum."
    cd ${dest_dir}
    md5sum -c `basename ${md5_file_path}`
    cd ${start_dir}

done

echo ""
echo ">> Successfully Exit pull.sh"
