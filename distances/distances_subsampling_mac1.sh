#!/bin/bash

vcf=$1
infix=subsampling/$(basename -s .vcf.gz $vcf)

for thin in ${infix}*.txt
do
	outfix=plink/$(basename -s .txt $thin)_mac1
	## PCA
	plink2 --vcf $vcf --keep $thin --out ${outfix}_pca --pca $(awk 'END{print NR-1}' $thin) --mac 1
	## FST
	plink2 --vcf $vcf --keep $thin --out ${outfix}_fstIID --fst family --pheno <(awk 'BEGIN{print "IID","family"} /deme/{$2 = $1; print}' $thin) --mac 1
	## DIST & MDS
	plink --vcf $vcf --keep <(sed 's/_/\t/' $thin) --out ${outfix}_dist --distance 1-ibs square0 --mac 1
	plink --vcf $vcf --keep <(sed 's/_/\t/' $thin) --out ${outfix}_mds --cluster --mds-plot 10 --mac 1
	## REL & COV MATRIX
	plink2 --vcf $vcf --keep $thin --out ${outfix}_rel --make-rel square0 --mac 1
	plink2 --vcf $vcf --keep $thin --out ${outfix}_cov --make-rel cov square0 --mac 1
done
