#!/bin/sh

#  picard.sh
#  Created by Kim Dill-McFarland on 10/13/21.

# Set vars
input="$1"
output="$2"
genome="$3"

if [ "$genome" == "Homo_sapiens.GRCh38" ]; then
    genome2="hg38"
elif [ "$genome" == "Mus_musculus.GRCm39" ]; then
    genome2="mm39"
else
    echo "Picard option currently only suports hg38 and mm39."
fi

# Download reference
if [ -f ref/PICARDref/refFlat.ensembl.txt ]; then
  echo "Picard reference already present in ref/PICARDref."
else
  mkdir -p ref/PICARDref
  sudo curl -O --output-dir ref/PICARDref \
    http://hgdownload.cse.ucsc.edu/goldenPath/${genome2}/database/refFlat.txt.gz
  gunzip ref/PICARDref/refFlat.txt.gz

  # Remove "chr" in chromosome name to match ensembl alignment using in STAR
  sed 's/chr//' ref/PICARDref/refFlat.txt > ref/PICARDref/refFlat.ensembl.txt
fi

# Run Picard
picard CollectRnaSeqMetrics \
    REF_FLAT=ref/PICARDref/refFlat.ensembl.txt \
    INPUT=${input} OUTPUT=${output} \
    ASSUME_SORTED=true STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=500 \
    QUIET=true VERBOSITY=ERROR
