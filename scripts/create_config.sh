#!/bin/bash

#################################
# Capture Samples to Config File
#################################

SampleList="data/*_R1*"

sudo echo "SampleList:" > result/config.yaml

spacer1=": "

# List samples & paired reads
for sample in $SampleList;
do
# Create R2 name
sample2=`echo "${sample/R1/R2}"`
# Create sample name
sample_name=`echo "$(basename $sample)" | grep -o '^.*_L' | sed 's/_L$//'`

# Add sample name to config
sudo echo "  " "$sample_name$spacer1" >> result/config.yaml
sudo echo "    sample: '"$sample_name"'" >> result/config.yaml
# Add fastq files to config
sudo echo "    R1: '"$sample"'" >> result/config.yaml
sudo echo "    R2: '"$sample2"'" >> result/config.yaml
done

#################################
# Add default param to Config File
#################################
# Auto detect cores
cores=$(eval nproc --all)
cores2=$(($cores-1))

# Add default param to config
sudo echo "

# Adapter removal
## Base pairs to trim from 5' end
trim5p: 10
## Removal of 3' adapter sequences? Default are Illumina Universal adapters
trimAdapt: True
adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# Genome alignment

## Species the format 'Homo_sapiens.GRCh38' or 'Mus_musculus.GRCm39'
genome: 'Homo_sapiens.GRCh38'
## Genome release number. Current as of 2022.04.14
release: '106'

# Alignment metrics
## Run Picard?
picard: True

# Other
threads: $cores2
" >> result/config.yaml

# Setup directory structure

sudo mkdir -p -m 777 'result/qc/1_fastqc_raw'
sudo mkdir -p -m 777 'result/qc/2_fastqc_trim'
sudo mkdir -p -m 777 'result/qc/3_flagstat'
sudo mkdir -p -m 777 'result/qc/4_picard'

sudo mkdir -p -m 777 'result/1_trim'
sudo mkdir -p -m 777 'result/2_bam'
sudo mkdir -p -m 777 'result/3_bam_filter'
sudo mkdir -p -m 777 'result/4_count'
sudo mkdir -p -m 777 'result/5_combined'

sudo mkdir -p -m 777 'ref'
sudo mkdir -p -m 777 'log'

