# mitogwas

A nextflow pipeline for running eQTL analysis of gene and RNA modifications encoded in mtDNA.

## Before using mitogwas:

This pipeline assumes that RNA sequencing data has been alinged to a reference genome with STAR (with --quantmode used), that RNAseQC has been run on each sample (with results in a subdirectory called RNAseQC) and that pileupAlleleExtractor_mito.pl (available in this repository) has been run on each sample, with an output name the same as the RNA name.

As an example, you should run the following for each RNA sequencing dataset (you can also include optional filtering steps after alignment):

```mkdir RNAseQC```

```STAR --runThreadN <num_threads> --genomeDir <path_to_genome_index> --readFilesIn <read1.fastq> [<read2.fastq>] --quantMode GeneCounts --outFileNamePrefix <output_prefix>```

```perl pileupAlleleExtractor_mito.pl --Bam output_prefix.Aligned.sortedByCoord.out.PP.UM.MT.bam --MinQ 23 --RefFasta <path_to_reference_fasta> --Out <output_prefix>```

```rnaseqc <path_to_gtf_file> output_prefix.Aligned.sortedByCoord.out.PP.UM.MT.bam RNAseQC```

## Pre-requististes:

Nextflow

Docker/Singularity

## Usage:

```nextflow run main.nf --rnaDir "RNAseqDIR" --bed "BEDFILE" --bim "BIMFILE" --fam "FAMFILE" --infoSheet "INFOSHEET" --gtfFile "GTFFILE" --dataName "NAME" -profile singularity/docker```

## Options:

--rnaDir "RNAseqDIR" : Full path to directory containing gene count files.  This directory should also contain a subdirectory called RNAseQC, which contains output of RNAseQC run on STAR alignment files.

--bed "BEDFILE" : Full path to plink format bed file, containing genetic variant data for samples to be used in the analysis.

--bim "BIMFILE" : Full path to plink format bim file, containing genetic variant data for samples to be used in the analysis.

--fam "FAMFILE" : Full path to plink format fam file, containing sample names (these should match those supplied in the inforsheet DNA column) to be used in the analysis.

--infoSheet "INFOSHEET": Full path to inforsheet file. This file is plain text with header "RNA DNA", followed by the names of any covariates to be used in the eQTL analysis. Below this there should be one line per sample to be included in the analysis, with the corresponding DNA (in fam file) and RNA (stub of STAR output) code for the sample, followed by values of any covarites to be used in the analysis.

--gtfFile "GTFFILE" : Full path to GTF file used for STAR aignment

--dataName "NAME" : Name for output

-profile singularity/docker : use either singularity or docker for the machine image

