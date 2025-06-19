use strict;
use Statistics::R;

#This program generates a normalised count matrices and QC plots for gene expression data in R

my $tail = $ARGV[0]; # e.g. Twins_Blood

my $R = Statistics::R->new();

#Counts distributions for batch:
my $make = "rawCountTable = read.table(\"batch_counts_all_tpm_$tail.txt\", header = TRUE, row.names = 1, check.names=FALSE)\ntest <- apply(rawCountTable, 2, function(x) { any(x > 0)})\nnew<-rawCountTable[test,]\nrawCountTable=new\npseudoCount = log10(rawCountTable + 1)\npdf(\"distributions_batch_counts_all_tpm_$tail.pdf\")\nhist(pseudoCount[,1],prob=TRUE,border=\"white\",main=\"RNAseq: TPM for all samples + log10(counts+1)\",xlab=\"TPM+log10(counts+1)\",ylim=c(0,4))\nfor (i in 1:ncol(pseudoCount)) { lines(density(pseudoCount[,i]),col=\"red\") }\ndev.off()";
$R->send($make);

#Count distributions to see outliers:
$make = "library(\"calibrate\")";
$R->send($make);
$make = "pdf(\"distributions_outliers_batch_counts_all_tpm_$tail.pdf\")\nhist(pseudoCount[,1],prob=TRUE,border=\"white\",main=\"RNAseq: TPM for all samples + log10(counts+1)\",xlab=\"TPM+log10(counts+1)\",ylim=c(0,200))\nfor (i in 1:ncol(pseudoCount)) { points(mean(pseudoCount[,i]),0.01,col=\"red\") }\nfor (i in 1:ncol(pseudoCount)) { textxy(mean(pseudoCount[,i]),(200-(i*0.075)),names(pseudoCount[i])) }\ndev.off()";
$R->send($make);

#PCA for Batch:
$make = "library(\"ggplot2\")";
$R->send($make);
$make = "pdf(\"pca_batch_counts_all_tpm_$tail.pdf\")\npca <- prcomp(t(pseudoCount),scale=T)\nsummary(pca)\nscores <- data.frame(pca\$x[,1:3])\nPoV <- pca\$sdev^2/sum(pca\$sdev^2)*100\na1=\"PC1: \"\na2=\"PC2: \"\na3=\"PC3: \"\nv1=round(PoV[1], digits=2)\nv2=round(PoV[2], digits=2)\nv3=round(PoV[3], digits=2)\nend=\"%\"\nt1=sprintf(\"\%s \%s\%s\",a1,v1,end)\nt2=sprintf(\"\%s \%s\%s\",a2,v2,end)\nt3=sprintf(\"\%s \%s\%s\",a3,v3,end)\np1<-qplot(x=PC1, y=PC2, data=scores, size=I(2),xlab = t1,ylab=t2)\np2<-qplot(x=PC2, y=PC3, data=scores, size=I(2),xlab = t2,ylab=t3)\np3<-qplot(x=PC1, y=PC3, data=scores, size=I(2),xlab = t1,ylab=t3)\nlibrary(gridExtra)\ngrid.arrange(p1, p2, p3, ncol=3)\ndev.off()";
$R->send($make);

#Open full TPM counts table, remove any inds with zero reads, remove genes that have zero reads in an individual, log transform, median norm, full quantile norm:
$make = "rawCountTable = read.table(\"counts_all_tpm_$tail.txt\", header = TRUE, row.names = 1, check.names=FALSE)\ntest <- apply(rawCountTable, 2, function(x) { any(x > 0)})\nnew<-rawCountTable[,test]\nrawCountTable=new\ntest <- apply(rawCountTable, 1, function(x) { all(x > 0)})\nnew<-rawCountTable[test,]\nrawCountTable=new\npseudoCount = log10(rawCountTable + 1)\nwrite.table(pseudoCount,file=\"log10_counts_all_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "pdf(\"distributions_log10_counts_all_tpm_$tail.pdf\")\nhist(pseudoCount[,1],prob=TRUE,border=\"white\",main=\"RNAseq: TPM for all samples + log10(counts+1)\",xlab=\"TPM+log10(counts+1)\",ylim=c(0,4))\nfor (i in 1:ncol(pseudoCount)) { lines(density(pseudoCount[,i]),col=\"red\") }\ndev.off()";
$R->send($make);

$make = "library(\"beadarray\")\nMedianpseudoCount=medianNormalise(pseudoCount, log=F)\nwrite.table(MedianpseudoCount,file=\"median_log10_counts_all_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "pdf(\"distributions_median_log10_counts_all_tpm_$tail.pdf\")\nhist(MedianpseudoCount[,1],prob=TRUE,border=\"white\",main=\"RNAseq: TPM for all samples + log10(counts+1) + Median\",xlab=\"TPM+log10(counts+1)\",ylim=c(0,4))\nfor (i in 1:ncol(MedianpseudoCount)) { lines(density(MedianpseudoCount[,i]),col=\"red\") }\ndev.off()";
$R->send($make);

$make = "library(\"preprocessCore\")\npseudonew=matrix(as.numeric(unlist(pseudoCount)),nrow=nrow(pseudoCount))\npseudoQ=normalize.quantiles(pseudonew)\nnames1=colnames(pseudoCount)\nnames2=rownames(pseudoCount)\ntry=as.data.frame(pseudoQ)\nnames(try) <- names1\nrow.names(try) <- names2\npseudoQ=try\nwrite.table(pseudoQ,file=\"quantile_log10_counts_all_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "pdf(\"distributions_quantile_log10_counts_all_tpm_$tail.pdf\")\nhist(pseudoQ[,1],prob=TRUE,border=\"white\",main=\"RNAseq: TPM for all samples + log10(counts+1) + Quantile\",xlab=\"TPM+log10(counts+1)\",ylim=c(0,4))\nfor (i in 1:ncol(pseudoQ)) { lines(density(pseudoQ[,i]),col=\"red\") }\ndev.off()";
$R->send($make);

#Open both mito (norm mito) counts table, remove any inds with zero reads, remove genes that have zero reads in an individual, log transform, median norm, full quantile norm:
$make = "rawCountTable = read.table(\"counts_mtmt_tpm_$tail.txt\", header = TRUE, row.names = 1, check.names=FALSE)\ntest <- apply(rawCountTable, 2, function(x) { any(x > 0)})\nnew<-rawCountTable[,test]\nrawCountTable=new\ntest <- apply(rawCountTable, 1, function(x) { all(x > 0)})\nnew<-rawCountTable[test,]\nrawCountTable=new\npseudoCount = log10(rawCountTable + 1)\nwrite.table(pseudoCount,file=\"log10_counts_mtmt_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "library(\"beadarray\")\nMedianpseudoCount=medianNormalise(pseudoCount, log=F)\nwrite.table(MedianpseudoCount,file=\"median_log10_counts_mtmt_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "library(\"preprocessCore\")\npseudonew=matrix(as.numeric(unlist(pseudoCount)),nrow=nrow(pseudoCount))\npseudoQ=normalize.quantiles(pseudonew)\nnames1=colnames(pseudoCount)\nnames2=rownames(pseudoCount)\ntry=as.data.frame(pseudoQ)\nnames(try) <- names1\nrow.names(try) <- names2\npseudoQ=try\nwrite.table(pseudoQ,file=\"quantile_log10_counts_mtmt_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

#Open both mito (norm nuc) counts table, remove any inds with zero reads, remove genes that have zero reads in an individual, log transform, median norm, full quantile norm:
$make = "rawCountTable = read.table(\"counts_mtnuc_tpm_$tail.txt\", header = TRUE, row.names = 1, check.names=FALSE)\ntest <- apply(rawCountTable, 2, function(x) { any(x > 0)})\nnew<-rawCountTable[,test]\nrawCountTable=new\ntest <- apply(rawCountTable, 1, function(x) { all(x > 0)})\nnew<-rawCountTable[test,]\nrawCountTable=new\npseudoCount = log10(rawCountTable + 1)\nwrite.table(pseudoCount,file=\"log10_counts_mtnuc_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "library(\"beadarray\")\nMedianpseudoCount=medianNormalise(pseudoCount, log=F)\nwrite.table(MedianpseudoCount,file=\"median_log10_counts_mtnuc_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

$make = "library(\"preprocessCore\")\npseudonew=matrix(as.numeric(unlist(pseudoCount)),nrow=nrow(pseudoCount))\npseudoQ=normalize.quantiles(pseudonew)\nnames1=colnames(pseudoCount)\nnames2=rownames(pseudoCount)\ntry=as.data.frame(pseudoQ)\nnames(try) <- names1\nrow.names(try) <- names2\npseudoQ=try\nwrite.table(pseudoQ,file=\"quantile_log10_counts_mtnuc_tpm_$tail.txt\",sep=\"\t\",col.names=NA, quote = FALSE)";
$R->send($make);

