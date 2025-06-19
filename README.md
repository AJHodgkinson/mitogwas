## mitogwas

A nextflow pipeline for running eQTL analysis of gene and RNA modifications encoded in mtDNA.

# Before using mitogwas:

This pipeline assumes that RNA sequencing data has been alinged to a reference genome with STAR (with --quantmode used), that RNAseQC has been run on each sample (with results in a subdirectory called RNAseQC) and that mitopileup.pl (available in this repository) has been run on each sample, with an output name the same as the RNA name.  As an example, if starting with paired data files: 

sample1.1.fastq.gz and sample1.2.fastq.gz

You should run:

mkdir RNAseQC

STAR

perl

Repeat this for every sample in your dataset, and then follow the instructions below to run mitogwas.

# Pre-requististes:

Nextflow

Docker/Sigularity

# Usage:

nextflow run main.nf --rnaDir "RNAseqDIR" --bed "BEDFILE" --bim "BIMFILE" --fam "FAMFILE" --infoSheet "INFOSHEET" --gtfFile "GTFFILE" --dataName "NAME" -profile singularity/docker

# Options:

--rnaDir "RNAseqDIR" : Full path to directory containing gene count files.  This directory should also contain a subdirectory called RNAseQC, which contains output of RNAseQC run on STAR alignment files.

--bed "BEDFILE" : Full path to plink format bed file, containing genetic variant data for samples to be used in the analysis.

--bim "BIMFILE" : Full path to plink format bim file, containing genetic variant data for samples to be used in the analysis.

--fam "FAMFILE" : Full path to plink format fam file, containing sample names (these should match those supplied in the inforsheet DNA column) to be used in the analysis.

--infoSheet "INFOSHEET": Full path to inforsheet file. This file is plain text with header "RNA DNA", followed by the names of any covariates to be used in the eQTL analysis. Below this there should be one line per sample to be included in the analysis, with the corresponding DNA (in fam file) and RNA (stub of STAR output) code for the sample, followed by values of any covarites to be used in the analysis.

--gtfFile "GTFFILE" : Full path to GTF file used for STAR aignment

--dataName "NAME" : Name for output

-profile singularity/docker : use either singularity or docker for the machine image

Pre-conditions:

Before running this pipeline, you should align your data to a reference genome using STAR. 
