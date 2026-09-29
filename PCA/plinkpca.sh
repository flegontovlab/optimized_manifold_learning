#!/bin/bash

vcf=$1
infix=subsampling/$(basename -s .vcf.gz $vcf)

for thin in ${infix}*.txt
do
    dims=$(awk 'END{print NR-1}' $thin)
    ## PCA with MAC=1
	outfix_mac=plink/$(basename -s .txt $thin)_mac1
	plink2 --vcf $vcf --keep $thin --out ${outfix_mac}_pca --pca ${dims} --mac 1

	## PCA with MAF=1%
	outfix_maf=plink/$(basename -s .txt $thin)_maf01
	plink2 --vcf $vcf --keep $thin --out ${outfix_maf}_pca --pca ${dims} --maf 0.01
done
