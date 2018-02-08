#!/bin/bash -eu

echo ">> Enter run_maf2maf.sh."

# Location of Java and Cromwell on your machine
java="/usr/bin/java"
cromwell_jar="/data/pipeline/v1/tools/cromwell/v0.21/cromwell-0.21.jar"

# Root Direcotry of git repository
root_dir=`cd ../../../../; pwd`
work_dir="${root_dir}/mutation_counts/examples/workflow/maf2maf/tmp"

# Source Copy
src_wdl_file="${root_dir}/mutation_counts/workflow/maf2maf.wdl"
src_inputs_json="${root_dir}/mutation_counts/examples/workflow/maf2maf/maf2maf_inputs.json"

# Working Copy
work_wdl_file=${work_dir}/`basename ${src_wdl_file}`
work_inputs_json=${work_dir}/`basename ${src_inputs_json}`

echo "## -------------------------------------------"
echo "## Parameters:"
echo "##   java: ${java}"
echo "##   cromwell_jar: ${java}"
echo "##   root_dir=${root_dir}"
echo "##   work_dir=${work_dir}"
echo "##   src_wdl_file=${src_wdl_file}"
echo "##   src_inputs_json=${src_inputs_json}"
echo "##   work_wdl_file=${work_wdl_file}"
echo "##   work_inputs_json=${work_inputs_json}"
echo "## -------------------------------------------"

if [ ! -f "${java}" ]; then
    echo ""
    echo ">> java '${java}' doesn't exist! Please set path to where Java is installed."
    exit 1
fi

if [ ! -f "${cromwell_jar}" ]; then
    echo ""
    echo ">> cromwell_jar '${cromwell_jar}' doesn't exist! Please set path to where Cromwell is installed."
    exit 1
fi


if [ -d "${work_dir}" ]; then
    echo ""
    echo ">> work_dir already exist. Removing it!"
    rm -v -r -f ${work_dir}
fi

echo ""
echo ">> Create workdir"
mkdir -v -p ${work_dir}

echo ""
echo ">> Create Working Copy of WDL and Inputs_Json files in work_dir"
cp ${src_wdl_file} ${work_wdl_file}
cp ${src_inputs_json} ${work_inputs_json}

echo ""
echo ">> Replace '<ROOT_DIR>' placeholder with actual value of root_dir."
# Note: Since '/' char appear in filepath, use the  '|' char as delimiter instead.
sed -i -e 's|<ROOT_DIR>|'"${root_dir}"'|g' ${work_inputs_json}

echo ""
echo ">> Change directory to work_dir."
cd ${work_dir}

echo ""
echo ">> Execute Cromwell."
${java} -jar ${cromwell_jar} run ${work_wdl_file} ${work_inputs_json}

echo ""
echo ">> Exit run_maf2maf.sh."
