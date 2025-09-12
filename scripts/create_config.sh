#!/bin/bash

#################################
# Capture Samples to Config File
#################################

SampleList="data/*_R1*"

echo "SampleList:" > result/config.yaml

spacer1=": "

# List samples & paired reads
for sample in $SampleList;
do
#Error is no files found
if [[ "$sample" == "data/*fastq.gz" || "$sample" == "data/*_R1*" ]]; then
echo "Error: No files found. Please check fastq naming requirements in the SEAsnake vignette."
break
fi

# Create R2 name
sample2=`echo "${sample/R1/R2}"`
# Create sample name
sample_name=`echo "$(basename $sample)" | grep -o '^.*_L[0-9][0-9][0-9]' | sed 's/_L[0-9][0-9][0-9]$//'`

# Add sample name to config
echo "  " "$sample_name$spacer1" >> result/config.yaml
echo "    sample: '"$sample_name"'" >> result/config.yaml
# Add fastq files to config
echo "    R1: '"$sample"'" >> result/config.yaml
echo "    R2: '"$sample2"'" >> result/config.yaml
done

#################################
# Add default param to Config File
#################################
# Auto detect cores
cores=$(eval nproc --all)
cores2=$(($cores-1))
#Fix is cores < 1
if [[ "$cores2" -lt 1 ]]; then
cores2=1
fi

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

## Species the format 'Homo_sapiens.GRCh38' or 'Mus_musculus.GRCm39'
genome: 'Homo_sapiens.GRCh38'
## Genome release number. Current as of 2025.06.05
release: '114'

# Alignment metrics
## Run Picard?
picard: True

# Other
threads: $cores2
" >> result/config.yaml

# Setup directory structure

mkdir -p -m 777 'result/qc/1_fastqc_raw'
mkdir -p -m 777 'result/qc/2_fastqc_trim'
mkdir -p -m 777 'result/qc/3_flagstat'
mkdir -p -m 777 'result/qc/4_picard'

mkdir -p -m 777 'result/1_trim'
mkdir -p -m 777 'result/2_bam'
mkdir -p -m 777 'result/3_bam_filter'
mkdir -p -m 777 'result/4_count'
mkdir -p -m 777 'result/5_combined'

mkdir -p -m 777 'ref'
mkdir -p -m 777 'log'

