#!/bin/bash

#################################
# Capture Samples to Config File
#################################


samples_list=(data/samples/*)


echo "samples:" > config/test_config.yaml
spacer=": "
indent="    "

# List samples
for sample in $samples_list;
do 
sample_name_simple="$(basename $sample)"
sample_name_simple2=`echo "$sample_name_simple" | cut -d'.' -f1`


sample_name_spaced="$indent$sample_name_simple2$spacer$sample_name_simple2"
echo $sample_name_spaced >> config/test_config.yaml
done

echo "extension:
    .html: .html
    _fastqc.zip: _fastqc.zip" >> config/test_config.yaml




