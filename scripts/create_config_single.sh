#!/bin/bash

#################################
# Capture Samples to Config File
#################################

SampleList="data/*fastq.gz"

sudo echo "SampleList:" > result/config.yaml

spacer1=": "

# List samples & paired reads
for sample in $SampleList;
do
# Create sample name
sample_name=`echo "$(basename $sample)" | sed 's/.fastq.gz$//'`

# Add sample name to config
sudo echo "  " "$sample_name$spacer1" >> result/config.yaml
sudo echo "    sample: '"$sample_name"'" >> result/config.yaml
# Add fastq files to config
sudo echo "    R1: '"$sample"'" >> result/config.yaml
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

## Species in one of the following formats:
### 'Homo_sapiens.GRCh38' for human
### 'Mus_musculus.GRCm39' for mouse
### an NCBI genome assembly like 'GCF_000195955.2_ASM19595v2' (Mtb)
genome: 'Homo_sapiens.GRCh38'

## Genome release number. Current as of 2024.12.06
## Ignored if using GCF genome
release: '113'

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

