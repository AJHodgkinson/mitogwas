#!/usr/bin/env nextflow

params.rnaDir = ''
params.bed = ''
params.bim = ''
params.fam = ''
params.infoSheet = ''
params.gtfFile = ''
params.dataName = ''

params.outDir = './results'

process qcFiles {

    input:
    path input

    output:
    path "failed_qc.txt"
    
    script:
    """
    perl ${baseDir}/bin/rnaseqc.pl $input >> failed_qc.txt
    """
}

process makeMatrix {

    publishDir "$params.outDir/", pattern: 'counts_all_tpm*', mode: 'copy'

    input:
    path input
    path failed_qc
    path infosheet
    val dataname
    path gtffile

    output:
    path "counts_all_raw_${dataname}.txt"
    path "gene_lengths.txt"
    path "counts_all_tpm_${dataname}.txt"
    path "batch_counts_all_tpm_${dataname}.txt"
    path "counts_mtnuc_tpm_${dataname}.txt"
    path "counts_mtmt_tpm_${dataname}.txt"
    path "meth_${dataname}.txt"

    script:
    """
    perl ${baseDir}/bin/htseq_create_matrix_nf.pl $input $failed_qc $infosheet counts_all_raw_${dataname}.txt
    perl ${baseDir}/bin/htseq_create_matrix_tpm_nf.pl $input $failed_qc $infosheet counts_all_tpm_${dataname}.txt counts_mtnuc_tpm_${dataname}.txt $gtffile
    perl ${baseDir}/bin/htseq_create_matrix_tpm_MT_nf.pl $input $failed_qc $infosheet counts_mtmt_tpm_${dataname}.txt $gtffile
    perl ${baseDir}/bin/methylation_matrix_nf.pl $input $failed_qc $infosheet meth_${dataname}.txt 

    """
}

process convertMatrix {

    publishDir "$params.outDir/", pattern: 'distributions_outliers_batch*', mode: 'copy'
    publishDir "$params.outDir/", pattern: 'pca_batch_counts_all_tpm_*', mode: 'copy'
    publishDir "$params.outDir/", pattern: 'distributions_log10_counts_all_tpm*', mode: 'copy'
    publishDir "$params.outDir/", pattern: 'distributions_median_log10_counts_all_tpm*', mode: 'copy'


    input:
    val dataname
    path counts_raw
    path gene_len
    path counts_all
    path batch_count
    path counts_mtnuc
    path count_mtmt
    path meth

    output:
    path "log10_counts_all_tpm_${dataname}.txt"
    path "median_log10_counts_all_tpm_${dataname}.txt"
    path "quantile_log10_counts_all_tpm_${dataname}.txt"
    path "log10_counts_mtmt_tpm_${dataname}.txt"
    path "median_log10_counts_mtmt_tpm_${dataname}.txt"
    path "quantile_log10_counts_mtmt_tpm_${dataname}.txt"
    path "log10_counts_mtnuc_tpm_${dataname}.txt"
    path "median_log10_counts_mtnuc_tpm_${dataname}.txt"
    path "quantile_log10_counts_mtnuc_tpm_${dataname}.txt"
    path "distributions_batch_counts_all_tpm_${dataname}.pdf"
    path "distributions_outliers_batch_counts_all_tpm_${dataname}.pdf"
    path "pca_batch_counts_all_tpm_${dataname}.pdf"
    path "distributions_log10_counts_all_tpm_${dataname}.pdf"
    path "distributions_median_log10_counts_all_tpm_${dataname}.pdf"
    path "distributions_quantile_log10_counts_all_tpm_${dataname}.pdf"

    script:
    """
    perl ${baseDir}/bin/convert_matrices_nf.pl ${dataname}

    """
}

process createPEER {
    input:
    val dataname
    path log10_all
    
    output:
    path "log10_counts_all_tpm_${dataname}.peer.csv"
    path "log10_counts_all_tpm_${dataname}.peer.header.txt"
    path "peer_nocovars_n50_log10_${dataname}"
    
    script:
    """
    perl ${baseDir}/bin/create_peer_exp_file_nf.pl $log10_all
    peertool -f log10_counts_all_tpm_${dataname}.peer.csv -n 50 -o peer_nocovars_n50_log10_${dataname}
    """
}

process createGPC {

    publishDir "$params.outDir/", pattern: '*maf5miss1hwe*', mode: 'copy'


    input:
    val dataname
    path bed
    path bim
    path fam
    path log10_all
    path infosheet

    output:
    path "${dataname}_maf5miss1hwe.bed"
    path "${dataname}_maf5miss1hwe.bim"
    path "${dataname}_maf5miss1hwe.fam"
    path "smart_${dataname}_merged_pruned12.out.evec"
    
    script:
    """
    awk 'NR==1 {for (i=1; i<=NF; i++) print \$i}' $log10_all > rna.txt
    grep -wf rna.txt $infosheet | cut -f1 > dna.txt 
    grep -wf dna.txt ${fam} | cut -d " " -f1,2 > keep.id.txt
    plink --bed $bed --bim $bim --fam $fam --keep keep.id.txt --maf 0.05 --geno 0.01 --make-bed --out ${dataname}_maf5miss1hwe --noweb
    awk 'NR==1 {print; next} { \$2 = \$1 } 1' ${dataname}_maf5miss1hwe.fam > ${dataname}_maf5miss1hwe.hold
    mv ${dataname}_maf5miss1hwe.hold ${dataname}_maf5miss1hwe.fam
    plink --bfile ${dataname}_maf5miss1hwe --indep-pairwise 100 5 0.3  --out ${dataname}_merged_prune --noweb
    plink --bfile ${dataname}_maf5miss1hwe  --extract ${dataname}_merged_prune.prune.in --recode12 --out ${dataname}_merged_pruned12 --noweb
    awk '{\$2=1 ; print ;}' ${dataname}_merged_pruned12.ped > ${dataname}_merged_hold.ped
    mv ${dataname}_merged_hold.ped ${dataname}_merged_pruned12.ped
    awk '{ if (length(\$2)>39) { print \$2 }}' ${dataname}_merged_pruned12.map > ${dataname}_merged_remove_long_indels.txt
    plink --file ${dataname}_merged_pruned12 --exclude ${dataname}_merged_remove_long_indels.txt --noweb --recode --out ${dataname}_merged_pruned12.short
    smartpca.perl -i ${dataname}_merged_pruned12.short.ped -a ${dataname}_merged_pruned12.short.map -b ${dataname}_merged_pruned12.short.ped -k 10 -o smart_${dataname}_merged_pruned12.out -p smart_${dataname}_merged_pruned12.par -e smart_${dataname}_merged_pruned12.eval -l smart_${dataname}_merged_pruned12.log -m 0	

    """

}


process makePLINK {

    publishDir "$params.outDir/", mode: 'copy'

    input:
    val dataname
    path meth_file
    path log_mtnuc
    path log_mtmt
    path med_mtnuc
    path med_mtmt
    path quant_mtnuc
    path quant_mtmt
    path infosheet
    path peer_path
    path peer_header
    path gpc_file

    output:
    path "plink_pheno_${dataname}.txt"
    path "plink_pheno_${dataname}_unmasked.txt"
    path "*pdf"

    script:
    """
    perl ${baseDir}/bin/make_plink_files_nf.pl $dataname $meth_file $log_mtnuc $log_mtmt $med_mtnuc $med_mtmt $quant_mtnuc $quant_mtmt $infosheet ${peer_path}/X.csv $peer_header $gpc_file
    """
}

process logexpASSOC {

    publishDir "$params.outDir/", pattern: '*assoc.linear.gz', mode: 'copy'
    publishDir "$params.outDir/", pattern: '*jpeg', mode: 'copy'

    input:
    val dataname
    path infosheet
    path bed
    path bim
    path fam
    path plink_pheno

    output:
    path "hold1.txt"
    path "*assoc.linear.gz"
    path "*jpeg"

    script:
    """
    store=\$(head -n 1 $infosheet | awk '{\$1=""; \$2=""; sub("  ", " "); print}' | cut -c 2- | sed 's/\s/,/g')    
    covars=\$(echo "\$store,GPC1,GPC2,GPC3,GPC4,GPC5,PEER1,PEER2,PEER3,PEER4,PEER5,PEER6,PEER7,PEER8,PEER9,PEER10")

    awk -v covar=\$covars 'NR==1 {
    for (i=1; i<=NF; i++) {
    	if ((\$i~/ENSG/)&&(\$i~/mtnuc/)&&(\$i~/log/)) {

	cmd = "plink --bfile ${dataname}_maf5miss1hwe --pheno $plink_pheno --pheno-name " \$i " --covar $plink_pheno --covar-name " covar " --linear hide-covar --no-const-covar --ci 0.95 --allow-no-sex --out ${dataname}." \$i
    	system(cmd)
	 }
      }
    }' $plink_pheno

    echo "hold" > hold1.txt
    echo "hello" > protect.jpeg
    for i in *linear; do stub=\$(echo \$i | sed 's/.assoc.linear//'); echo -e "library(\\"qqman\\")\nresults = read.table(\\"\$i\\", header=T)\njpeg(\\"\${stub}.jpeg\\",quality=100)\nqq(results\\\$P)\ndev.off()" > \${stub}.R; done
    for i in *R; do stub=\$(echo \$i | sed 's/\\.R//'); val=\$(awk 'NR==2 { print \$12 }' \$stub.assoc.linear); [[ \$val =~ [0-9] ]] && Rscript \$i; done
    for i in *linear; do pigz \$i; done

    """
    
}

process medexpASSOC {

    publishDir "$params.outDir/", pattern: '*assoc.linear.gz', mode: 'copy'
    publishDir "$params.outDir/", pattern: '*jpeg', mode: 'copy'

    input:
    val dataname
    path infosheet
    path bed
    path bim
    path fam
    path plink_pheno
    path hold1

    output:
    path "hold2.txt"
    path "*assoc.linear.gz"
    path "*jpeg"

    script:
    """
    store=\$(head -n 1 $infosheet | awk '{\$1=""; \$2=""; sub("  ", " "); print}' | cut -c 2- | sed 's/\s/,/g')    
    covars=\$(echo "\$store,GPC1,GPC2,GPC3,GPC4,GPC5,PEER1,PEER2,PEER3,PEER4,PEER5,PEER6,PEER7,PEER8,PEER9,PEER10")

    awk -v covar=\$covars 'NR==1 {
    for (i=1; i<=NF; i++) {
    	if ((\$i~/ENSG/)&&(\$i~/mtnuc/)&&(\$i~/med/)) {

	cmd = "plink --bfile ${dataname}_maf5miss1hwe --pheno $plink_pheno --pheno-name " \$i " --covar $plink_pheno --covar-name " covar " --linear hide-covar --no-const-covar --ci 0.95 --allow-no-sex --out ${dataname}." \$i
    	system(cmd)
	 }
      }
    }' $plink_pheno

    echo "hold" > hold2.txt
    echo "hello" > protect.jpeg
    for i in *linear; do stub=\$(echo \$i | sed 's/.assoc.linear//'); echo -e "library(\\"qqman\\")\nresults = read.table(\\"\$i\\", header=T)\njpeg(\\"\${stub}.jpeg\\",quality=100)\nqq(results\\\$P)\ndev.off()" > \${stub}.R; done
    for i in *R; do stub=\$(echo \$i | sed 's/\\.R//'); val=\$(awk 'NR==2 { print \$12 }' \$stub.assoc.linear); [[ \$val =~ [0-9] ]] && Rscript \$i; done
    for i in *linear; do pigz \$i; done

    """
    
}

process methASSOC {

    publishDir "$params.outDir/", pattern: '*assoc.linear.gz', mode: 'copy'
    publishDir "$params.outDir/", pattern: '*jpeg', mode: 'copy'

    input:
    val dataname
    path infosheet
    path bed
    path bim
    path fam
    path plink_pheno
    path hold2

    output:
    path "*assoc.linear.gz"
    path "*jpeg"

    script:
    """
    store=\$(head -n 1 $infosheet | awk '{\$1=""; \$2=""; sub("  ", " "); print}' | cut -c 2- | sed 's/\s/,/g')    
    covars=\$(echo "\$store,GPC1,GPC2,GPC3,GPC4,GPC5,PEER1,PEER2,PEER3,PEER4,PEER5,PEER6,PEER7,PEER8,PEER9,PEER10")

    awk -v covar=\$covars 'NR==1 {
    for (i=1; i<=NF; i++) {
    	if (\$i~/_full/) {

	cmd = "plink --bfile ${dataname}_maf5miss1hwe --pheno $plink_pheno --pheno-name " \$i " --covar $plink_pheno --covar-name " covar " --linear hide-covar --no-const-covar --ci 0.95 --allow-no-sex --out ${dataname}." \$i
    	system(cmd)
	 }
      }
    }' $plink_pheno

    echo "hello" > protect.jpeg
    for i in *linear; do stub=\$(echo \$i | sed 's/.assoc.linear//'); echo -e "library(\\"qqman\\")\nresults = read.table(\\"\$i\\", header=T)\njpeg(\\"\${stub}.jpeg\\",quality=100)\nqq(results\\\$P)\ndev.off()" > \${stub}.R; done
    for i in *R; do stub=\$(echo \$i | sed 's/\\.R//'); val=\$(awk 'NR==2 { print \$12 }' \$stub.assoc.linear); [[ \$val =~ [0-9] ]] && Rscript \$i; done
    for i in *linear; do pigz \$i; done

    """
    
}



workflow {   
    qc_ch = qcFiles(params.rnaDir)
    matrix_ch = makeMatrix(params.rnaDir, qc_ch, params.infoSheet, params.dataName, params.gtfFile)
    convert_ch = convertMatrix(params.dataName, matrix_ch)
    peer_ch = createPEER(params.dataName,convert_ch.getAt(0))
    gpc_ch = createGPC(params.dataName,params.bed,params.bim,params.fam,convert_ch.getAt(0),params.infoSheet)
    plink_ch = makePLINK(params.dataName,matrix_ch.getAt(6),convert_ch.getAt(6),convert_ch.getAt(3),convert_ch.getAt(7),convert_ch.getAt(4),convert_ch.getAt(8),convert_ch.getAt(5),params.infoSheet,peer_ch.getAt(2),peer_ch.getAt(1),gpc_ch.getAt(3))
    log_assoc_ch = logexpASSOC(params.dataName, params.infoSheet, gpc_ch.getAt(0), gpc_ch.getAt(1), gpc_ch.getAt(2), plink_ch.getAt(0))
    med_assoc_ch = medexpASSOC(params.dataName, params.infoSheet, gpc_ch.getAt(0), gpc_ch.getAt(1), gpc_ch.getAt(2), plink_ch.getAt(0), log_assoc_ch.getAt(0))
    meth_assoc_ch = methASSOC(params.dataName, params.infoSheet, gpc_ch.getAt(0), gpc_ch.getAt(1), gpc_ch.getAt(2), plink_ch.getAt(0), med_assoc_ch.getAt(0))
}












