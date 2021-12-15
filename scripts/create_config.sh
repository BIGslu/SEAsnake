#!/bin/bash

# Setup directory structure
mkdir -p 'result/qc/1_fastqc_raw'
mkdir -p 'result/qc/2_fastqc_trim'
mkdir -p 'result/qc/3_flagstat'
mkdir -p 'result/qc/4_picard'

mkdir -p 'result/1_trim'
mkdir -p 'result/2_bam'
mkdir -p 'result/3_bam_filter'
mkdir -p 'result/4_count'
mkdir -p 'result/5_combined'

mkdir -p 'ref'
mkdir -p 'log'

#################################
# Capture Samples to Config File
#################################

SampleList="data/*R1*"


echo "SampleList:" > result/config.yaml

spacer1=": "

# List samples & paired reads
for sample in $SampleList;
do
# Create R2 name
sample2=`echo "${sample/R1/R2}"`
# Create sample name
sample_name=`echo "$(basename $sample)" | grep -o '^.*_L' | sed 's/_L$//'`

# Add sample name to config
echo "  " "$sample_name$spacer1" >> result/config.yaml
echo "    sample: '"$sample_name"'" >> result/config.yaml
# Add fastq files to config
echo "    R1: '"$sample"'" >> result/config.yaml
echo "    R2: '"$sample2"'" >> result/config.yaml
done

# Auto detect cores
cores=$(eval nproc --all)
cores2=$(($cores-1))

# Add default param to config
echo "
# Adapter removal
## Base pairs to trim from 5' end
trim5p: 10
## Removal of 3' adapter sequences? Default are Illumina Universal adapters
trimAdapt: True
adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# Genome alignment
## Human genome release number
release: '104'

# Alignment metrics
## Run Picard?
picard: True

# Other
threads: $cores2
" >> result/config.yaml
