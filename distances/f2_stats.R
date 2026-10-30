library(admixtools)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
samples <- args[1] # subsampling/*.txt

prefix <- paste0("fstats/", basename(samples)) %>% str_replace("gens.+", "gens")
outfix <- paste0("fstats/", basename(samples), "_f2")

print(prefix)

pops <- read_csv(samples, col_names = FALSE)
pops <- pops %>% transmute(pop = str_remove(X1, "_."), ind = X1)

extract_f2(
  pref = prefix,
  outdir = outfix,
  inds = pops$ind,
  blgsize = 100000000,
  fst = FALSE,
  afprod = FALSE,
  #minmaf = 0.01,
  adjust_pseudohaploid = FALSE,
  overwrite = TRUE,
  verbose = TRUE,
  n_cores = 32
)
fstats <- f2_from_precomp(outfix)

f2s <- f2(fstats, sure = TRUE, verbose = TRUE)
write_tsv(f2s, paste0(outfix, ".tsv"))
