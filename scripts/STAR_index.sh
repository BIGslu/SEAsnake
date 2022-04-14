#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
species="$2"
threads="$3"

if [ $species == "human" ]; then
    full_species="homo_sapiens"
    genome="Homo_sapiens.GRCh38"
fi

if[ $species == "mouse" ]; then
    full_species="mus_musculus"
    genome="Mus_musculus.GRCm39"
fi

# Setup directories
mkdir -p ref/release${release}/STARref
mkdir -p ref/release${release}/STARindex

# Download gtf if not present
gtf=ref/release${release}/STARref/${genome}.${release}.gtf

if [ ! -e "$gtf" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/gtf/${full_species}/${genome}.${release}.gtf.gz
    yes y | gunzip ref/release${release}/STARref/*gtf.gz
else
    echo "Genome GTF already exists. No new file downloaded."
fi

# Download fasta if not present
fasta=ref/release${release}/STARref/${genome}.dna.primary_assembly.fa

if [ ! -e "$fasta" ]; then
    sudo curl -O --output-dir ref/release${release}/STARref ftp://ftp.ensembl.org/pub/release-${release}/fasta/${full_species}/dna/${genome}.dna.primary_assembly.fa.gz
    yes y | gunzip ref/release${release}/STARref/*fa.gz
else
    echo "Genome fasta already exists. No new file downloaded."
fi

# Make index if not present
SA=ref/release${release}/STARindex/SA

if [ ! -e "$SA" ]; then
    outDir=ref/release${release}/STARindex
    
    STAR --runMode genomeGenerate --genomeDir ${outDir} --genomeFastaFiles ref/release${release}/STARref/${genome}.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/${genome}.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
else
    echo "Genome index already exists. No new index created."
fi

# Paste genome name to config
echo "
genome: '"$genome"'
">> result/config.yaml