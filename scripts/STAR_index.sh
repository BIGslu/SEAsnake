#!/bin/sh

#  STAR_index.sh
#  Created by Kim Dill-McFarland on 10/14/21.
#

# Set vars
release="$1"
genome="$2"
threads="$3"

#Setup for human and mouse
if [[ "$genome" == *"Homo_sapiens"* || "$genome" == *"Mus_musculus"* ]]; then
  species=(${genome//./ })
  species="${species%% *}"
  species=`echo "$species" | awk '{print tolower($0)}'`
  #dir
  out=ref/release${release}
  ref=${out}/STARref
  index=${out}/STARindex
  SA=${out}/STARindex/SA
  #url
  fasta_url=ftp://ftp.ensembl.org/pub/release-${release}/fasta/${species}/dna/${genome}.dna.primary_assembly.fa.gz
  gtf_url=ftp://ftp.ensembl.org/pub/release-${release}/gtf/${species}/${genome}.${release}.gtf.gz
  #files
  gtf=${out}/STARref/${genome}.${release}.gtf
  fasta=${out}/STARref/${genome}.dna.primary_assembly.fa
#Setup for GCF genomes
elif [[ "$genome" == "GCF_"* ]]; then
  species=$genome
  #dir
  out=ref/${genome}
  ref=${out}/STARref
  index=${out}/STARindex
  SA=${out}/STARindex/SA
  #Parse GCF number
  numbers="${genome#*_}" 
  numbers="${numbers%%.*}" 
  v1="${numbers:0:3}"
  v2="${numbers:3:3}"
  v3="${numbers:6:3}" 
  #url
  fasta_url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${v1}/${v2}/${v3}/${genome}/${genome}_genomic.fna.gz
  gtf_url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${v1}/${v2}/${v3}/${genome}/${genome}_genomic.gtf.gz
  #files
  gtf=${out}/STARref/${genome}_genomic.gtf
  fasta=${out}/STARref/${genome}_genomic.fna
else
  echo "Error: Genome input not recognized."
  exit 1
fi

# Setup directories
if [ ! -d "$ref" ]; then
    mkdir -p ${ref}
fi

if [ ! -d "$index" ]; then
    mkdir -p ${index}
fi

# Download gtf if not present
if [ ! -e "$gtf" ]; then
    sudo curl -s -O --output-dir ${out}/STARref ${gtf_url}
    yes y | gunzip ${out}/STARref/*gtf.gz
else
  echo "Genome GTF already exists. No new file downloaded."
fi

# Download fasta if not present
if [ ! -e "$fasta" ]; then
    sudo curl -s -O --output-dir ${out}/STARref ${fasta_url}
    yes y | gunzip ${out}/STARref/*.gz
else
  echo "Genome fasta already exists. No new file downloaded."
fi

# Make index if not present
if [ ! -e "$SA" ]; then
    STAR --runMode genomeGenerate --genomeDir ${index} --genomeFastaFiles {$fasta} --sjdbGTFfile {$gtf} --sjdbOverhang 99 --runThreadN ${threads}
else
    echo "Genome index already exists. No new index created."
fi
