#!/bin/bash

#################################
# Capture Samples to Config File
#################################


SampleList="data/*R1*"


echo "SampleList:" > config/test_config.yaml

spacer1=": "
spacer2=""
indent="  "

# List samples & paired reads
for sample in $SampleList;
do 
sample_name1_simple="$(basename $sample)"
sample_name=`echo "$sample_name1_simple" | grep -o '^.*_L' | sed 's/_L$//'`

sample_name1_simple2=`echo "$sample_name1_simple" | cut -d'.' -f1`

sample_name2_simple2=`echo "${sample_name1_simple2/R1/R2}"`

# sample_name_spaced="$indent$sample_name$spacer1$sample_name1_simple2$spacer2$sample_name2_simple2"
echo "  " "$sample_name$spacer1" >> config/test_config.yaml
echo "    sample: '"$sample_name"'" >> config/test_config.yaml
# echo $sample_name_spaced >> config/test_config.yaml
echo "    R1: 'data/"$sample_name1_simple2"'" >> config/test_config.yaml
echo "    R2: 'data/"$sample_name2_simple2"'" >> config/test_config.yaml
done

# Build the rest of the config file
echo "
trim5p: 10
trimAdapt: True
adapter1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
release: 104
threads: 6
" >> config/test_config.yaml




