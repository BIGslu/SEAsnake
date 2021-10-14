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

# Download gtf if not present
gtf=ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf

if [ ! -e "$gtf" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.gtf.gz
    yes y | gunzip ref/release${release}/STARref/*gtf.gz
else
    echo "Genome GTF already exists. No new file downloaded."
fi

# Download fasta if not present
fasta=ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

if [ ! -e "$fasta" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    yes y | gunzip ref/release${release}/STARref/*fa.gz
else
    echo "Genome fasta already exists. No new file downloaded."
fi

# Make index if empty
if [ -z "$(ls -A ref/release${release}/STARindex)" ]; then
outDir=ref/release${release}/STARindex
echo $outDir
    #STAR --runMode genomeGenerate --genomeDir  --genomeFastaFiles ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
else
   echo "Files exist in ref/release#/STARindex. No new index created."
fi

