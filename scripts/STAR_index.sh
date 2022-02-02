#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
threads="$2"

ref=ref/release${release}/STARref
index=ref/release${release}/STARindex
gtf=ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf
fasta=ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SA=ref/release${release}/STARindex/SA

# Setup directories
if [ ! -d "$ref" ]; then
    mkdir -p ref/release${release}/STARref
fi

if [ ! -d "$index" ]; then
    mkdir -p ref/release${release}/STARindex
fi

# Download gtf if not present
if [ ! -e "$gtf" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.gtf.gz
    yes y | gunzip ref/release${release}/STARref/*gtf.gz
else
    echo "Genome GTF already exists. No new file downloaded."
fi

# Download fasta if not present
if [ ! -e "$fasta" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    yes y | gunzip ref/release${release}/STARref/*fa.gz
else
    echo "Genome fasta already exists. No new file downloaded."
fi

# Make index if not present
if [ ! -e "$SA" ]; then  
    STAR --runMode genomeGenerate --genomeDir ${index} --genomeFastaFiles ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
else
    echo "Genome index already exists. No new index created."
fi
