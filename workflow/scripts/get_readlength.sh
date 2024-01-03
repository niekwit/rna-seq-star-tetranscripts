#!/usr/bin/env bash

# read length - 1 (for STAR --sjdbOverhang) from MultiQC file
# if multiple read lengths are present, take the longest one
# if average read length is a float, round down to nearest integer
# check if read length -1 is an integer

set +e

awk '{print $5}' ${snakemake_input[0]} | sed 1d | sort -nr | uniq | head -1 | awk '{print $1-1}' | cut -d. -f1 | tee ${snakemake_output[0]} | awk '{if ($1 ~ /^[0-9]+$/) exit 0; else exit 1}' 2> ${snakemake_log[0]}

exitcode=$?

if [ $exitcode -eq 1 ]
then
    echo "ERROR: Read length -1 is not a number!"
    echo "Please check results/qc/multiqc_data/multiqc_general_stats.txt..."
    exit 1
else 
    exit 0
fi