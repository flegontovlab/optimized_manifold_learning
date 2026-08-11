#reticulate::use_python("/home/jena/miniconda3/envs/phate/bin/python", required = TRUE)
#library(phateR)
#library(umap)
library(uwot)
library(tidyverse)
library(magrittr)
library(future)
library(future.apply)

## test phate is working
# ph = phate(tree.data.small$data)
# str(ph)

#### PARAMETERS ####

# options(future.globals.maxSize = 1000*1024^2)
future::plan('multicore')
plan(list(tweak(multicore, workers = 48)))

## UMAP
dist_method <- c("euclidean", "manhattan", "cosine", "correlation")
## PHATE
# dist_method <- c("euclidean", "cityblock", "cosine")

## tidyr::expand_grid sorts by first column
## base::expand.grid sorts by last column
combs = c(seq(3,20,1), seq(22,40,2), seq(45,100,5), seq(110,199,10), 199) %>% expand.grid(.,.) %>% set_colnames(c("PCs", "knn")) #%>% glimpse

## annotation of populations
anno = read_csv("hexagon_landscape_colours.csv")

## ARGS = filename of samplelist
args = commandArgs(trailingOnly=TRUE)
prefix = args[1] # subsampling/longdistgfs_331demes_200mb_01_300pairs_10gens_clusters_0gen_200ind_thin1.txt

pcafiles = c(
	paste0("input/", prefix, "_center_evec.tsv"),
	paste0("input/", prefix, "_drift_evec.tsv"),
	paste0("input/", prefix, "_pca_alldims.eigenvec")
)

## samples are the same for the whole script
iidfile = paste0("input/", prefix, ".txt") #%>% str_remove("_maf01")
samples = read_tsv(iidfile, col_names = "IID") %>% arrange(IID) %>% mutate(FID = str_remove(IID, "_."))
## geo dist
coords = inner_join(samples, anno[,1:4], by = "FID")
dist_geo = coords %>% select(3:5) %>% dist %>% as.vector

## importing genetic distances
## vectors are faster than dataframes
dist_ibs <- read_tsv(paste0("input/", prefix,"_dist.mdist"), col_names=F) %>% as.dist %>% as.vector
dist_rel <- read_tsv(paste0("input/", prefix,"_rel.rel"), col_names=F) %>% as.dist %>% as.vector
dist_cov <- read_tsv(paste0("input/", prefix,"_cov.rel"), col_names=F) %>% as.dist %>% as.vector
dist_mds <- read.table(paste0("input/", prefix,"_mds.mds"), h = T)[,4:6] %>% dist %>% as.vector
dist_fst <- read_tsv(paste0("input/", prefix,"_fstIID.fst.summary"))[[3]]
dist_f2 <- read_tsv(paste0("input/", prefix,"_f2.tsv"))[[3]]
## fix negative values
dist_fst[dist_fst < 0] <- 0
dist_f2[dist_f2 < 0] <- 0

#### LOOPS ####

## loop over pca setups
for (pcafile in pcafiles){
	pca = read_tsv(pcafile) %>%
		select(where(is.numeric))
	dist_3PCs = pca %>% select(PC1,PC2,PC3) %>% dist %>% as.vector

	## get replicate info
	replicate = pcafile %>%
		basename %>% tools::file_path_sans_ext() %>%
		str_match(".*0mb_([0-9]{2})_.*_([[:lower:]]+)_.*_(thin.)_(.*)") %>%
		extract(2:5) %>%
		paste(collapse = "_")

	results = combs
	results$replicate = replicate

	## caclulating R2 vs. 3D PC space
	## needs to be done once per pca
	results$R2_geogr_Pearson_PC123 <- cor(dist_3PCs, dist_geo, method="pearson")^2
	results$R2_geogr_Spearman_PC123 <- cor(dist_3PCs, dist_geo, method="spearman")^2
	results$R2_IBS_def_Pearson_PC123 <- cor(dist_3PCs, dist_ibs, method="pearson")^2
	results$R2_IBS_def_Spearman_PC123 <- cor(dist_3PCs, dist_ibs, method="spearman")^2
	results$R2_rel_Pearson_PC123 <- cor(dist_3PCs, dist_rel, method="pearson")^2
	results$R2_rel_Spearman_PC123 <- cor(dist_3PCs, dist_rel, method="spearman")^2
	results$R2_cov_Pearson_PC123 <- cor(dist_3PCs, dist_cov, method="pearson")^2
	results$R2_cov_Spearman_PC123 <- cor(dist_3PCs, dist_cov, method="spearman")^2
	results$R2_mds_Pearson_PC123 <- cor(dist_3PCs, dist_mds, method="pearson")^2
	results$R2_mds_Spearman_PC123 <- cor(dist_3PCs, dist_mds, method="spearman")^2
	results$R2_FST_Pearson_PC123 <- cor(dist_3PCs, dist_fst, method="pearson")^2
	results$R2_FST_Spearman_PC123 <- cor(dist_3PCs, dist_fst, method="spearman")^2
	results$R2_f2_Pearson_PC123 <- cor(dist_3PCs, dist_f2, method="pearson")^2
	results$R2_f2_Spearman_PC123 <- cor(dist_3PCs, dist_f2, method="spearman")^2

	## loop over dist metrics
	for (d in dist_method){
		## running UMAP optimization, p1:
		rep <- c(1:nrow(combs))
		outtab = future_lapply(rep, function(rep){
			latent <- uwot::umap(pca[,1:combs$PCs[rep]], n_components=3, metric=d, n_neighbors=combs$knn[rep], min_dist=0.1, n_threads=1, verbose=F)
			dist_umap <- latent %>% dist %>% as.vector
			data.frame(
				dist_method = d,
				PCs = combs$PCs[rep],
				knn = combs$knn[rep],
				R2_geogr_Pearson_UMAP123 = cor(dist_umap, dist_geo, method = "pearson")^2,
				R2_geogr_Spearman_UMAP123 = cor(dist_umap, dist_geo, method = "spearman")^2,
				R2_IBS_def_Pearson_UMAP123 = cor(dist_umap, dist_ibs, method = "pearson")^2,
				R2_IBS_def_Spearman_UMAP123 = cor(dist_umap, dist_ibs, method = "spearman")^2,
				R2_rel_Pearson_UMAP123 = cor(dist_umap, dist_rel, method = "pearson")^2,
				R2_rel_Spearman_UMAP123 = cor(dist_umap, dist_rel, method = "spearman")^2,
				R2_cov_Pearson_UMAP123 = cor(dist_umap, dist_cov, method = "pearson")^2,
				R2_cov_Spearman_UMAP123 = cor(dist_umap, dist_cov, method = "spearman")^2,
				R2_mds_Pearson_UMAP123 = cor(dist_umap, dist_mds, method = "pearson")^2,
				R2_mds_Spearman_UMAP123 = cor(dist_umap, dist_mds, method = "spearman")^2,
				R2_FST_Pearson_UMAP123 = cor(dist_umap, dist_fst, method = "pearson")^2,
				R2_FST_Spearman_UMAP123 = cor(dist_umap, dist_fst, method = "spearman")^2,
				R2_f2_Pearson_UMAP123 = cor(dist_umap, dist_f2, method = "pearson")^2,
				R2_f2_Spearman_UMAP123 = cor(dist_umap, dist_f2, method = "spearman")^2
				)
		}) %>% bind_rows %>% inner_join(results, by = c("PCs", "knn"))
		write_tsv(outtab, file = paste0("umap/UMAP_full_optimization_", replicate, "_", d, "_3D.tsv"))

		## running UMAP optimization, p2
		outtab = future_lapply(rep, function(rep){
			latent <- uwot::umap(pca[,1:combs$PCs[rep]], n_components=3, metric=d, n_neighbors=combs$knn[rep], min_dist=0.5, n_threads=1, verbose=F)
			dist_umap <- latent %>% dist %>% as.vector
			data.frame(
				dist_method = d,
				PCs = combs$PCs[rep],
				knn = combs$knn[rep],
				R2_geogr_Pearson_UMAP123 = cor(dist_umap, dist_geo, method = "pearson")^2,
				R2_geogr_Spearman_UMAP123 = cor(dist_umap, dist_geo, method = "spearman")^2,
				R2_IBS_def_Pearson_UMAP123 = cor(dist_umap, dist_ibs, method = "pearson")^2,
				R2_IBS_def_Spearman_UMAP123 = cor(dist_umap, dist_ibs, method = "spearman")^2,
				R2_rel_Pearson_UMAP123 = cor(dist_umap, dist_rel, method = "pearson")^2,
				R2_rel_Spearman_UMAP123 = cor(dist_umap, dist_rel, method = "spearman")^2,
				R2_cov_Pearson_UMAP123 = cor(dist_umap, dist_cov, method = "pearson")^2,
				R2_cov_Spearman_UMAP123 = cor(dist_umap, dist_cov, method = "spearman")^2,
				R2_mds_Pearson_UMAP123 = cor(dist_umap, dist_mds, method = "pearson")^2,
				R2_mds_Spearman_UMAP123 = cor(dist_umap, dist_mds, method = "spearman")^2,
				R2_FST_Pearson_UMAP123 = cor(dist_umap, dist_fst, method = "pearson")^2,
				R2_FST_Spearman_UMAP123 = cor(dist_umap, dist_fst, method = "spearman")^2,
				R2_f2_Pearson_UMAP123 = cor(dist_umap, dist_f2, method = "pearson")^2,
				R2_f2_Spearman_UMAP123 = cor(dist_umap, dist_f2, method = "spearman")^2
				)
		}) %>% bind_rows %>% inner_join(results, by = c("PCs", "knn"))
		write_tsv(outtab, file = paste0("umap/UMAP_full_optimization_", replicate, "_", d, "_3D_mindist05.tsv"))

		## running UMAP optimization, p3
		outtab = future_lapply(rep, function(rep){
			latent <- uwot::umap(pca[,1:combs$PCs[rep]], n_components=3, metric=d, n_neighbors=combs$knn[rep], min_dist=0.1, dens_scale=0.5, n_trees=100, n_epochs=1000, n_threads=1, verbose=F)
			dist_umap <- latent %>% dist %>% as.vector
			data.frame(
				dist_method = d,
				PCs = combs$PCs[rep],
				knn = combs$knn[rep],
				R2_geogr_Pearson_UMAP123 = cor(dist_umap, dist_geo, method = "pearson")^2,
				R2_geogr_Spearman_UMAP123 = cor(dist_umap, dist_geo, method = "spearman")^2,
				R2_IBS_def_Pearson_UMAP123 = cor(dist_umap, dist_ibs, method = "pearson")^2,
				R2_IBS_def_Spearman_UMAP123 = cor(dist_umap, dist_ibs, method = "spearman")^2,
				R2_rel_Pearson_UMAP123 = cor(dist_umap, dist_rel, method = "pearson")^2,
				R2_rel_Spearman_UMAP123 = cor(dist_umap, dist_rel, method = "spearman")^2,
				R2_cov_Pearson_UMAP123 = cor(dist_umap, dist_cov, method = "pearson")^2,
				R2_cov_Spearman_UMAP123 = cor(dist_umap, dist_cov, method = "spearman")^2,
				R2_mds_Pearson_UMAP123 = cor(dist_umap, dist_mds, method = "pearson")^2,
				R2_mds_Spearman_UMAP123 = cor(dist_umap, dist_mds, method = "spearman")^2,
				R2_FST_Pearson_UMAP123 = cor(dist_umap, dist_fst, method = "pearson")^2,
				R2_FST_Spearman_UMAP123 = cor(dist_umap, dist_fst, method = "spearman")^2,
				R2_f2_Pearson_UMAP123 = cor(dist_umap, dist_f2, method = "pearson")^2,
				R2_f2_Spearman_UMAP123 = cor(dist_umap, dist_f2, method = "spearman")^2
				)
		}) %>% bind_rows %>% inner_join(results, by = c("PCs", "knn"))
		write_tsv(outtab, file = paste0("umap/UMAP_full_optimization_", replicate, "_", d, "_3D_densscale05_ntrees100_nepochs1000.tsv"))

		## running UMAP optimization, p4
		outtab = future_lapply(rep, function(rep){
			latent <- uwot::umap(pca[,1:combs$PCs[rep]], n_components=3, metric=d, n_neighbors=combs$knn[rep], min_dist=0.1, dens_scale=1, n_trees=100, n_epochs=1000, n_threads=1, verbose=F)
			dist_umap <- latent %>% dist %>% as.vector
			data.frame(
				dist_method = d,
				PCs = combs$PCs[rep],
				knn = combs$knn[rep],
				R2_geogr_Pearson_UMAP123 = cor(dist_umap, dist_geo, method = "pearson")^2,
				R2_geogr_Spearman_UMAP123 = cor(dist_umap, dist_geo, method = "spearman")^2,
				R2_IBS_def_Pearson_UMAP123 = cor(dist_umap, dist_ibs, method = "pearson")^2,
				R2_IBS_def_Spearman_UMAP123 = cor(dist_umap, dist_ibs, method = "spearman")^2,
				R2_rel_Pearson_UMAP123 = cor(dist_umap, dist_rel, method = "pearson")^2,
				R2_rel_Spearman_UMAP123 = cor(dist_umap, dist_rel, method = "spearman")^2,
				R2_cov_Pearson_UMAP123 = cor(dist_umap, dist_cov, method = "pearson")^2,
				R2_cov_Spearman_UMAP123 = cor(dist_umap, dist_cov, method = "spearman")^2,
				R2_mds_Pearson_UMAP123 = cor(dist_umap, dist_mds, method = "pearson")^2,
				R2_mds_Spearman_UMAP123 = cor(dist_umap, dist_mds, method = "spearman")^2,
				R2_FST_Pearson_UMAP123 = cor(dist_umap, dist_fst, method = "pearson")^2,
				R2_FST_Spearman_UMAP123 = cor(dist_umap, dist_fst, method = "spearman")^2,
				R2_f2_Pearson_UMAP123 = cor(dist_umap, dist_f2, method = "pearson")^2,
				R2_f2_Spearman_UMAP123 = cor(dist_umap, dist_f2, method = "spearman")^2
				)
		}) %>% bind_rows %>% inner_join(results, by = c("PCs", "knn"))
		write_tsv(outtab, file = paste0("umap/UMAP_full_optimization_", replicate, "_", d, "_3D_densscale1_ntrees100_nepochs1000.tsv"))
	}
}
