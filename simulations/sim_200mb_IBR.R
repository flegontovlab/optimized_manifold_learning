library(tidyverse)
library(magrittr)
library(slendr)
init_env()

## args
args = commandArgs(trailingOnly = T)
simrep = args[1] # filledcircle01

outprefix = paste(simrep, "331demes_200mb", sep = "_")

outpath = "./sim/"
dir.create(outpath)
outprefix = paste0(outpath, outprefix)
print(outprefix)

## landscape and neighbours
landscape <- read_csv("landscape.csv")
neighbour_pairs <- read_tsv("neighbour_pairs.tsv")

popsize = 1000

## adm. prop. intervals
preLGMmin = 0.006
preLGMmax = 1
LGMmin = 0.0002
LGMmax = 0.25
postLGMmin = 0.003
postLGMmax = 1

n1 <- round(rnorm(10000, mean=0, sd=0.3), digits=3)
n2 <- round(rnorm(10000, mean=0, sd=0.1), digits=4)
n3 <- round(rnorm(10000, mean=0, sd=0.3), digits=3)
df1 <- as.data.frame(sample(n1[n1>preLGMmin & n1<preLGMmax], size = nrow(neighbour_pairs))) # 326
df2 <- as.data.frame(sample(n2[n2>LGMmin & n2<LGMmax], size = nrow(neighbour_pairs))) # 326
df3 <- as.data.frame(sample(n3[n3>postLGMmin & n3<postLGMmax], size = nrow(neighbour_pairs))) # 326
colnames(df1) <- c("gene_flows")
colnames(df2) <- c("gene_flows")
colnames(df3) <- c("gene_flows")

i = 20
j = 10
preLGMdur = 1400
LGMdur = 350
postLGMdur = 650

preLGMst = j+23*i+1
preLGMend = j+23*i+preLGMdur
LGMst = preLGMend+1
LGMend = preLGMend+LGMdur
postLGMst = LGMend+1
postLGMend = j+23*i+preLGMdur+LGMdur+6*i+postLGMdur

p_anc <- population("p_anc", N=1000, time=1)

## non-spatial populations
## pop list
pops <- vector("list", length = nrow(landscape))
names(pops) <- landscape$FID

for (i in 1:length(pops)) {
    pops[[i]] = population(names(pops[i]), N = 1000, time = 470, map = F, parent = p_anc)
}

## concatenate all pops
poplist <- c(list(p_anc), pops)

## define gene flows
## empty lists per epoch
gf_preLGM1 <- vector("list", length = nrow(neighbour_pairs))
gf_preLGM2 <- vector("list", length = nrow(neighbour_pairs))
gf_preLGM3 <- vector("list", length = nrow(neighbour_pairs))
gf_preLGM4 <- vector("list", length = nrow(neighbour_pairs))
gf_preLGM5 <- vector("list", length = nrow(neighbour_pairs))
gf_LGM1 <- vector("list", length = nrow(neighbour_pairs))
gf_LGM2 <- vector("list", length = nrow(neighbour_pairs))
gf_LGM3 <- vector("list", length = nrow(neighbour_pairs))
gf_LGM4 <- vector("list", length = nrow(neighbour_pairs))
gf_LGM5 <- vector("list", length = nrow(neighbour_pairs))
gf_postLGM1 <- vector("list", length = nrow(neighbour_pairs))
gf_postLGM2 <- vector("list", length = nrow(neighbour_pairs))
gf_postLGM3 <- vector("list", length = nrow(neighbour_pairs))
gf_postLGM4 <- vector("list", length = nrow(neighbour_pairs))
gf_postLGM5 <- vector("list", length = nrow(neighbour_pairs))

## fill the lists in a loop
for (i in 1:nrow(neighbour_pairs)) {
    ## preLGM
    gf_preLGM1[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df1[i,1], preLGMst, preLGMst+1*preLGMdur/5-1, overlap=F)
    gf_preLGM2[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df1[i,1], preLGMst+1*preLGMdur/5, preLGMst+2*preLGMdur/5-1, overlap=F)
    gf_preLGM3[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df1[i,1], preLGMst+2*preLGMdur/5, preLGMst+3*preLGMdur/5-1, overlap=F)
    gf_preLGM4[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df1[i,1], preLGMst+3*preLGMdur/5, preLGMst+4*preLGMdur/5-1, overlap=F)
    gf_preLGM5[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df1[i,1], preLGMst+4*preLGMdur/5, preLGMst+preLGMdur-1, overlap=F)
    ## LGM
    gf_LGM1[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df2[i,1], LGMst, LGMst+1*LGMdur/5-1, overlap=F)
    gf_LGM2[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df2[i,1], LGMst+1*LGMdur/5, LGMst+2*LGMdur/5-1, overlap=F)
    gf_LGM3[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df2[i,1], LGMst+2*LGMdur/5, LGMst+3*LGMdur/5-1, overlap=F)
    gf_LGM4[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df2[i,1], LGMst+3*LGMdur/5, LGMst+4*LGMdur/5-1, overlap=F)
    gf_LGM5[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df2[i,1], LGMst+4*LGMdur/5, LGMst+LGMdur-1, overlap=F)
    ## postLGM
    gf_postLGM1[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df3[i,1], postLGMst, postLGMst+1*(postLGMend-LGMend)/5-1, overlap=F)
    gf_postLGM2[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df3[i,1], postLGMst+1*(postLGMend-LGMend)/5, postLGMst+2*(postLGMend-LGMend)/5-1, overlap=F)
    gf_postLGM3[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df3[i,1], postLGMst+2*(postLGMend-LGMend)/5, postLGMst+3*(postLGMend-LGMend)/5-1, overlap=F)
    gf_postLGM4[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df3[i,1], postLGMst+3*(postLGMend-LGMend)/5, postLGMst+4*(postLGMend-LGMend)/5-1, overlap=F)
    gf_postLGM5[[i]] <- gene_flow(pops[[neighbour_pairs$pop1[i]]], pops[[neighbour_pairs$pop2[i]]], df3[i,1], postLGMst+4*(postLGMend-LGMend)/5, postLGMst+(postLGMend-LGMend)-1, overlap=F)
}

## concatenate all gf lists
gf <- c(
    gf_preLGM1, gf_preLGM2, gf_preLGM3, gf_preLGM4, gf_preLGM5,
    gf_LGM1, gf_LGM2, gf_LGM3, gf_LGM4, gf_LGM5,
    gf_postLGM1, gf_postLGM2, gf_postLGM3, gf_postLGM4, gf_postLGM5
)

## compile the model and define additional parameters
model <- compile_model(
  populations = poplist,
  gene_flow = gf, generation_time = 1, simulation_length = postLGMend,
  serialize = T, overwrite = T, force = T, path = paste0(outprefix, "_scripts")
  )

## schedule sampling
## the "schedule_sampling" function cannot be used programmatically,
## but its output is just a table, which *can* be made programmatically!
sampling_plan <- data.frame(
  time = c(postLGMend-300, postLGMend),
  pop = rep(landscape$FID, each = 2),
  n = as.integer(3),
  x = NA, y = NA, x_orig = NA, y_orig = NA
)

######## SIMULATE ########
## this code is for slendr 0.8 or 0.9
ts <- msprime(model, sequence_length=200e6, recombination_rate=1e-8, samples=sampling_plan, verbose=T, output=paste0(outprefix, ".ts"))

## process the results
ts_coalesced(ts)
ts <- ts_mutate(ts, mutation_rate=1.25e-8)

write_tsv(ts_samples(ts), file="filledcircle_multifurcation_samples.txt")

#snps <- ts_eigenstrat(ts, prefix=outprefix, chrom = "1")

ts_vcf(ts, path = paste0(outprefix, ".vcf.gz"), chrom = "1")

#sessionInfo()
