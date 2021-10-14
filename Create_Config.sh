#!/bin/bash

#################################
# Capture Samples to Config File
#################################


samples_list="data/samples/*R1*"


echo "sample_list:" > config/test_config.yaml

spacer1=": "
spacer2=" "
indent="    "

# List samples & paired reads
for sample in $samples_list;
do 
sample_name1_simple="$(basename $sample)"
sample_name=`echo "$sample_name1_simple" | grep -o '^.*_L' | sed 's/_L$//'`

sample_name1_simple2=`echo "$sample_name1_simple" | cut -d'.' -f1`

sample_name2_simple2=`echo "${sample_name1_simple2/R1/R2}"`

# sample_name_spaced="$indent$sample_name$spacer1$sample_name1_simple2$spacer2$sample_name2_simple2"
echo "  " "$sample_name$spacer1" >> config/test_config.yaml
echo "    sample: " $sample_name >> config/test_config.yaml
# echo $sample_name_spaced >> config/test_config.yaml
echo "    R1: " $sample_name1_simple2 >> config/test_config.yaml
echo "    R2: " $sample_name2_simple2 >> config/test_config.yaml
done

# Build the rest of the config file
echo "
extension:
    .html: .html
    _fastqc.zip: _fastqc.zip
	
trim5p:
    10
trimadapt:
    True
adapter1:
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2:
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
release:
    104
" >> config/test_config.yaml




