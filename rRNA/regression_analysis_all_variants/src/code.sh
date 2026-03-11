#!/bin/bash

#run with --priority ---destnation rDNA\ Variations:regression_analysis --instance-type mem2_ssd1_v2_x4
#regression_type is logistic or linear
#dx run -ichunk=0 -iregression_type=logistic -y regression_analysis --priority low --destination rDNA\ Variations:regression_analysis --instance-type mem1_hdd1_v2_x16
# Main script to download references, prepare environment, and run single_run.sh in parallel

# Reference files and directories
script="regression_analysis_on_chunk.py"
mkdir -p /home/dnanexus/output_folder

#coding_path="/home/dnanexus/mnt/rDNA Variations/coding19.tsv"
#field="/home/dnanexus/mnt/rDNA Variations/field.txt"
#covariates_path="/home/dnanexus/mnt/rDNA Variations/covariates.extra.csv"
#variants_path="/home/dnanexus/mnt/rDNA Variations/merged_rdna_variant_frequencies.unrelated_WB.csv"
#disease_path="/home/dnanexus/mnt/rDNA Variations/icd10_and_cancers.all.csv"


dx download rDNA\ Variations:coding19.tsv
dx download rDNA\ Variations:field.txt
dx download rDNA\ Variations:covariates.extra.csv
dx download rDNA\ Variations:merged_rdna_variant_frequencies.unrelated_WB.csv
dx download rDNA\ Variations:merged_rdna_count_frequencies.all.csv
dx download rDNA\ Variations:icd10_and_cancers.all.csv
dx download rDNA\ Variations:phenotypes.csv

set -e -x
apt-get update
apt-get install -y --fix-missing samtools bowtie2 python3 git wget tar

chmod u+x /usr/bin/dxfuse
#dx-mount-all-inputs
#ls -l in

mkdir /home/dnanexus/mnt

#echo "Checking if mount works and we can see the input files"
dxfuse /home/dnanexus/mnt "rDNA Variations"

#chmod u+x ./single_run.sh
log_file="command_logs.txt"
touch $log_file

past_run_output_path="/home/dnanexus/mnt/rDNA Variations/regression_analysis_all_variants/"

output_files=()
###ran with chunk size = 25 for control_type = have_no_diseases, chunk size = 5 for control_type = any, chunk size = 5 for linear

# Determine chunk size based on control_type
if [ "$control_type" == "have_no_diseases" ]; then
  chunk_size=25
elif [ "$control_type" == "any" ]; then
  chunk_size=5
else
  chunk_size=5  # Default value if control_type is unspecified or different
fi

start=$((chunk * chunk_size))
end=$((start + chunk_size - 1))

output_chunk="disease_${regression_type}_converge_reg.${chunk}.${control_type}.with_genetic_covariates.unrelated_WB.csv"

if [ -e "${past_run_output_path}/${output_chunk}" ]; then
  echo "disease_${regression_type}_converge_reg.${chunk}.${control_type}.with_genetic_covariates.unrelated_WB.csv File exists." | tee -a $log_file
else
  echo "Running custom Python script for ${regression_type} ${chunk}" | tee -a $log_file
  python3 "$script" "${regression_type}" "${start}" "${end}" "${output_chunk}" "${control_type}"
  echo "Created ${output_chunk}" | tee -a $log_file
  file_id=$(dx upload "$output_chunk" --brief)
  dx-jobutil-add-output output_files "$file_id" --class=array:file
fi
#
#file_id=$(dx upload "${log_file}" --brief)
#dx-jobutil-add-output output_files "$file_id" --class=array:file

