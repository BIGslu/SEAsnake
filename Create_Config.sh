#!/bin/bash

#################################
# Capture Samples to Config File
#################################


samples_list="data/samples/*R1*"

#echo $samples_list

echo "samples:" > config/test_config.yaml
spacer1=": "
spacer2=" "
indent="    "

# List samples
for sample in $samples_list;
do 
sample_name1_simple="$(basename $sample)"
sample_name=`echo "$sample_name1_simple" | grep -o '^.*_L' | sed 's/_L$//'`

sample_name1_simple2=`echo "$sample_name1_simple" | cut -d'.' -f1`

sample_name2_simple2=`echo "${sample_name1_simple2/R1/R2}"`

sample_name_spaced="$indent$sample_name$spacer1$sample_name1_simple2$spacer2$sample_name2_simple2"
echo $sample_name
echo $sample_name_spaced >> config/test_config.yaml
done

echo "extension:
    .html: .html
    _fastqc.zip: _fastqc.zip" >> config/test_config.yaml




