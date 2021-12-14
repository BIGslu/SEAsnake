#!/bin/sh

#  picard.sh
#  Created by Kim Dill-McFarland on 10/13/21.

# Set vars
input="$1"
output="$2"

# Download reference
if [ ref/PICARDref/refFlat.ensembl.txt ]; then
  echo "Picard reference already present in ref/PICARDref."
else
  mkdir -p ref/PICARDref
  sudo curl -O --output-dir ref/PICARDref \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
  gunzip ref/PICARDref/refFlat.txt.gz

  # Remove "chr" in chromosome name to match ensembl alignment using in STAR
  sed 's/chr//' ref/PICARDref/refFlat.txt > ref/PICARDref/refFlat.ensembl.txt

fi

# Run Picard
 java -jar ~/apps/anaconda/share/picard-2.26.2-0/picard.jar \
    CollectRnaSeqMetrics \
    REF_FLAT=ref/PICARDref/refFlat.ensembl.txt \
    INPUT=${input} OUTPUT=${output} \
    ASSUME_SORTED=true STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=500 \
    QUIET=true VERBOSITY=ERROR
