#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
output="$2"
threads="$3"

# Setup directories
mkdir -p ref/release${release}/STARref
mkdir -p ref/release${release}/STARindex

# Make ftp path names
gtf=ftp://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.gtf.gz
fasta=ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download reference
sudo curl -O --output-dir ref/release${release}/STARref ${gtf}
sudo curl -O --output-dir ref/release${release}/STARref ${fasta}
gunzip ref/release${release}/STARref/*

# Index genome
STAR --runMode genomeGenerate --genomeDir ${output} --genomeFastaFiles ref/release${release}/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ref/release${release}/STARref/Homo_sapiens.GRCh38.${release}.gtf --sjdbOverhang 99 --runThreadN ${threads}
