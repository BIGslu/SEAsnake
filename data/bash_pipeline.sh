#!/bin/bash

##### Bulk paired-end RNAseq data cleaning ######

########################################
## Setup
########################################

read1="data/test_S1_L005_R1_001.fastq.gz"
read2="data/test_S1_L005_R2_001.fastq.gz"
name="test"
threads=6
trim5p=10
adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
release=104

# Create directory structure
mkdir -p results/fastqc
mkdir -p results/fastq_trim
mkdir -p results/bam
mkdir -p results/metrics
mkdir -p results/counts
mkdir -p ref/

########################################
## Quality assessment 1
########################################

fastqc $read1 -o results/fastqc/ -t $threads
fastqc $read2 -o results/fastqc/ -t $threads

########################################
## Adapter removal
########################################

AdapterRemoval --file1 $read1 --file2 $read2 \
    --basename results/fastq_trim/$name --gzip \
    --trim5p $trim5p --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --adapter1 $adapter1 --adapter2 $adapter2 \
    --threads $threads

########################################
## Quality assessment 2
########################################

fastqc results/fastq_trim/$name.pair1.truncated.gz -o results/fastqc/ -t $threads
fastqc results/fastq_trim/$name.pair2.truncated.gz -o results/fastqc/ -t $threads

########################################
## Alignment reference
########################################

#Set file structure
mkdir -p ref/release$release/STARref

#Download reference
sudo curl -O --output-dir ref/release$release/STARref \
    ftp://ftp.ensembl.org/pub/release-$release/gtf/homo_sapiens/Homo_sapiens.GRCh38.$release.gtf.gz
    
sudo curl -O --output-dir ref/release$release/STARref \
    ftp://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip ref/release$release/STARref/*

#Index genome
STAR --runMode genomeGenerate \
     --genomeDir ref/release$release/STARindex \
     --genomeFastaFiles \
       ref/release$release/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile ref/release$release/STARref/Homo_sapiens.GRCh38.$release.gtf \
     --sjdbOverhang 99 \
     --runThreadN $threads

#Move log file to within reference directory
mv Log.out ref/release$number/STARindex/

########################################
## Alignment
########################################
STAR --genomeDir ref/release$release/STARindex \
         --readFilesIn fastq_trim/$name.pair1.truncated.gz fastq_trim/$name.pair2.truncated.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix results/bam/$name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $threads \
         --runRNGseed 8756

########################################
## Quality filter alignment
########################################

samtools view results/bam/$name.bam \
    -h -f 3 -F 1284 -q 30 -@ $threads > results/bam/"$name"_filter.bam

########################################
## flagstat
########################################

samtools flagstat -@ $threads results/bam/$name.bam > results/metrics/"$name"_flagstat.tsv
samtools flagstat -@ $threads results/bam/"$name"_filter.bam > results/metrics/"$name"_filter_flagstat.tsv

########################################
## Picard
########################################
#Set file structure
mkdir -p ref/PICARDref

#Download reference
sudo curl -O --output-dir ref/PICARDref \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
    
gunzip ref/PICARDref/refFlat.txt.gz

#Run Picard
java -jar ~/apps/anaconda/share/picard-2.26.2-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=ref/PICARDref/refFlat.ensembl.txt \
        INPUT=results/bam/"$name"_filter.bam \
        OUTPUT=results/metrics/"$name"_picard.tsv \
        ASSUME_SORTED=true STRAND_SPECIFICITY=NONE MINIMUM_LENGTH=500 \
        QUIET=true VERBOSITY=ERROR

########################################
## Count reads in genes
########################################

featureCounts -T $threads -g gene_id -t exon -p \
  -a ref/release$release/STARref/Homo_sapiens.GRCh38.$release.gtf \
  -o results/counts/"$name"_featurecounts.tsv \
  results/bam/*filter.bam

################# END ##################