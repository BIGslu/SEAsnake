# SEAsnake
RNA-seq pipeline

This pipeline includes quality assessment and filtering, alignment, and count table generation. Specifically,

1. Quality assess sequences with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Remove adapters and filter low quality sequences with [AdapterRemoval](https://adapterremoval.readthedocs.io/en/latest/)
3. Align to reference genome with [STAR](https://github.com/alexdobin/STAR)
4. Quality filter alignments with [samtools](http://www.htslib.org/) `view`
5. Quality assess alignments with [samtools](http://www.htslib.org/) `flagstat` and/or [Picard](https://broadinstitute.github.io/picard/) `CollectRnaSeqMetrics`
6. Count reads in genes with [Subread](http://subread.sourceforge.net/) `featureCounts`

See our [tutorial](https://github.com/BIGslu/tutorials/blob/main/RNAseq/1.Hawn_RNAseq_fastq.to.counts.pdf) for further details.
