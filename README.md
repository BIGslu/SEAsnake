# SEAsnake

[![DOI](https://zenodo.org/badge/411778345.svg)](https://zenodo.org/badge/latestdoi/411778345)

Citation: Dill-McFarland KA, Benson B, Segnitx RM. 2024. SEAsnake: a snakemake pipeline for processing RNA-seq fastq to counts. v1.1. Zenodo. doi: 10.5281/zenodo.11646755

## RNA-seq pipeline

This pipeline includes quality assessment and filtering, alignment, and count table generation. Specifically,

1. Quality assess sequences with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Remove adapters and filter low quality sequences with [AdapterRemoval](https://adapterremoval.readthedocs.io/en/latest/)
3. Align to reference genome with [STAR](https://github.com/alexdobin/STAR)
4. Quality filter alignments with [samtools](http://www.htslib.org/) `view`
5. Quality assess alignments with [samtools](http://www.htslib.org/) `flagstat` and/or [Picard](https://broadinstitute.github.io/picard/) `CollectRnaSeqMetrics`
6. Count reads in genes with [Subread](http://subread.sourceforge.net/) `featureCounts`

See our [tutorial](https://bigslu.github.io/SEAsnake/vignette/SEAsnake_vignette.html) for further details.
