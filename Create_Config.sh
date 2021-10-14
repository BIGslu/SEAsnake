#!/bin/bash

#################################
# Capture Samples to Config File
#################################


SampleList="data/*R1*"


echo "SampleList:" > config/test_config.yaml

spacer1=": "

# List samples & paired reads
for sample in $SampleList;
do
# Create R2 name
sample2=`echo "${sample/R1/R2}"`
# Create sample name
sample_name=`echo "$(basename $sample)" | grep -o '^.*_L' | sed 's/_L$//'`

# Add sample name to config
echo "  " "$sample_name$spacer1" >> config/test_config.yaml
echo "    sample: '"$sample_name"'" >> config/test_config.yaml
# Add fastq files to config
echo "    R1: '"$sample"'" >> config/test_config.yaml
echo "    R2: '"$sample2"'" >> config/test_config.yaml
done

# Add default param to config
echo "
trim5p: 10
trimAdapt: True
adapter1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
release: 104
picard: True
threads: 6
" >> config/test_config.yaml




