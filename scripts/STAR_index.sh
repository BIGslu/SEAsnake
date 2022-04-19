#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
genome="$2"
threads="$3"
species=(${genome//./ })
species=`echo "$species" | awk '{print tolower($0)}'`

ref=ref/release${release}/STARref
index=ref/release${release}/STARindex
gtf=ref/release${release}/STARref/${genome}.${release}.gtf
fasta=ref/release${release}/STARref/${genome}.dna.primary_assembly.fa
SA=ref/release${release}/STARindex/SA
    
# Setup directories
if [ ! -d "$ref" ]; then
    mkdir -p ${ref}
fi

if [ ! -d "$index" ]; then
    mkdir -p ${index}
fi

# Download gtf if not present
if [ ! -e "$gtf" ]; then
    sudo curl -s -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/gtf/${species}/${genome}.${release}.gtf.gz
    yes y | gunzip ref/release${release}/STARref/*gtf.gz
else
    echo "Genome GTF already exists. No new file downloaded."
fi

# Download fasta if not present
if [ ! -e "$fasta" ]; then
    sudo curl -s -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/fasta/${species}/dna/${genome}.dna.primary_assembly.fa.gz
    yes y | gunzip ref/release${release}/STARref/*fa.gz
else
    echo "Genome fasta already exists. No new file downloaded."
fi

# Make index if not present
if [ ! -e "$SA" ]; then
    STAR --runMode genomeGenerate --genomeDir ${index} --genomeFastaFiles ref/release${release}/STARref/${genome}.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/${genome}.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
else
    echo "Genome index already exists. No new index created."
fi
