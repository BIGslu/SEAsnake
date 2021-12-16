#!/bin/bash

#################################
# Run FastQC on raw fastq files
#################################

# Set vars
input="$1"

# Auto detect cores
cores=$(eval nproc --all)
cores2=$(($cores-1))

# FastQC
fastqc -t ${cores2} -o result/qc/1_fastqc_raw ${input}
