library(tidyverse)
library(magrittr)
library(smartsnp)

args <- commandArgs(trailingOnly = TRUE)
popfile <- args[1]
scaling <- args[2] # "drift" or "center"

datafile <- sub(".txt", ".traw", popfile)
outprefix <- sub(".txt", "", popfile) %>% basename %>% paste0("pca/", .)
print(outprefix)

poplist <- read_tsv(popfile, col_names = "IID") %>%
  arrange(IID) %>% # traw files were sorted by plink, my poplists were not
  mutate(FID = str_remove(IID, "_."))

pca <- smart_pca(snp_data = datafile,
                 sample_group = poplist$IID,
                 pc_axes = nrow(poplist),
                 scaling = scaling)

evec <- as.data.frame(pca$pca.sample_coordinates)
eval <- data.frame(eval = pca$pca.eigenvalues[1,])
eval$PC <- 1:nrow(poplist)
#loadings <- as.data.frame(pca$pca.snp_loadings)

write_tsv(evec, paste0(outprefix, "_", scaling, "_evec.tsv"))
write_tsv(eval, paste0(outprefix, "_", scaling, "_eval.tsv"))
#write_tsv(loadings, paste0(outprefix, "_", scaling, "_loadings.tsv"))
