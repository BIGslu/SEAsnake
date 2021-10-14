#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
threads="$2"

# Setup directories
mkdir -p ref/release${release}/STARref
mkdir -p ref/release${release}/STARindex

# Download reference if empty
if [ -z "$(ls -A ref/release${release}/STARref)" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.gtf.gz
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    # Unzip
    yes y | gunzip ref/release${release}/STARref/*
else
   echo "Files exist in ref/release#/STARref. No new files downloaded."
fi

# Make index if empty
if [ -z "$(ls -A ref/release${release}/STARindex)" ]; then
outDir=ref/release${release}/STARindex
echo $outDir
    #STAR --runMode genomeGenerate --genomeDir  --genomeFastaFiles ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
else
   echo "Files exist in ref/release#/STARindex. No new index created."
fi

