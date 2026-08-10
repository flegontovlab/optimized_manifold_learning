#!/bin/bash

## sim rep
for vcf in sim/*200mb*.vcf.gz # {01..10}
do
	outname=$(basename -s .vcf.gz $vcf)
	## thin rep
	for j in {1..5}
	do
		## clustered sampling
		## 0gen epoch
		sort -R demes.txt | head -5 | grep -f - neighbours.tsv | tr "\t" "\n" | sort -u | grep -f - samples_0gen.txt | sort -R | head -50 > subsampling/${outname}_clusters_0gen_50ind_thin${j}.txt
		sort -R demes.txt | head -15 | grep -f - neighbours.tsv | tr "\t" "\n" | sort -u | grep -f - samples_0gen.txt | sort -R | head -200 > subsampling/${outname}_clusters_0gen_200ind_thin${j}.txt
		## 300gen epoch
		sort -R demes.txt | head -5 | grep -f - neighbours.tsv | tr "\t" "\n" | sort -u | grep -f - samples_300gen.txt | sort -R | head -50 > subsampling/${outname}_clusters_300gen_50ind_thin${j}.txt
		sort -R demes.txt | head -15 | grep -f - neighbours.tsv | tr "\t" "\n" | sort -u | grep -f - samples_300gen.txt | sort -R | head -200 > subsampling/${outname}_clusters_300gen_200ind_thin${j}.txt

		## random sampling
		## 0gen epoch
		sort -R samples_0gen.txt | head -50 > subsampling/${outname}_random_0gen_50ind_thin${j}.txt
		sort -R samples_0gen.txt | head -200 > subsampling/${outname}_random_0gen_200ind_thin${j}.txt
		## 300gen epoch
		sort -R samples_300gen.txt | head -50 > subsampling/${outname}_random_300gen_50ind_thin${j}.txt
		sort -R samples_300gen.txt | head -200 > subsampling/${outname}_random_300gen_200ind_thin${j}.txt
	done
done
