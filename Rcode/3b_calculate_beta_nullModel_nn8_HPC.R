###################################################################################################################
## calculate expected phylogenetic and functional beta-diversity given taxonomic beta-diversity within regions of 3*3 grid-cells 
## by randomize species in specific regions in the phylogeny or trait dendrogram  
###################################################################################################################

rm(list = ls())

# Set working directory
path2wd <- "/work/wubing/global_tree_beta_2022"
setwd(path2wd)

# load packages
needed_libs <- c("betapart", "raster", "spdep", "picante", "ape", "phytools", "parallel", "abind")

for(x in needed_libs){
    if(!require(x, character.only=TRUE)){
			if(!require(x,character.only=TRUE, lib.loc = "/gpfs0/home/wubing/R/library")){
				install.packages(x, repos = "http://cran.us.r-project.org", lib.loc = "/gpfs0/home/wubing/R/library", dependencies = TRUE)
				require(x, character.only = TRUE, lib.loc = "/gpfs0/home/wubing/R/library")					
			}
  }
}

source("Rcode/functions/get_beta.R")

load("data/tree_pam/tree_pam6_final.RDATA")
load("data/traits_phylogeny/trait_phylogeny_final.RDATA")

# community matrix
comdata <- tree_pam6$Pre

# set parallel computing
t1 <- Sys.time()
print("n_cores")
print(detectCores())
no_cores <- 25   # detectCores() the number of cores
cl <- makeCluster(no_cores, type="PSOCK") # Initiate cluster
clusterExport(cl, varlist = c("comdata", "tree", "trait.tree", "get_beta"))
clusterEvalQ(cl, c(library(betapart,), library(spdep, lib.loc = "/gpfs0/home/wubing/R/library"), 
			library(picante), library(ape), library(phytools), library(abind)))

## phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=comdata, tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=5, d2=300, ncell_resample=5, null="local"), silent=TRUE)
phylo.sor_nn8_null <- parLapply(cl,1:200, get_beta_tempo)
phylo.sor_nn8_null <- sapply(phylo.sor_nn8_null, as.matrix, simplify="array")

save(phylo.sor_nn8_null, file="intermediate_results/random_beta_phylo_nn8_raw.RDATA")

## functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=comdata, tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=5, d2=300, ncell_resample=5, null="local"), silent=TRUE)
func.sor_nn8_null <- parLapply(cl,1:200, get_beta_tempo)
func.sor_nn8_null <- sapply(func.sor_nn8_null, as.matrix, simplify="array")

save(phylo.sor_nn8_null, func.sor_nn8_null, file="intermediate_results/random_beta_phylo_func_nn8_raw.RDATA")

stopCluster(cl)
t2 <- Sys.time()
t2 - t1

## calculate mean and sd of returned measures
## phylogenetic beta
# add maximum and minnimum alpha and their differences and percentage of nestedness
phylo.sor_nn8_null <- abind(phylo.sor_nn8_null, 
                        "alpha.max" = phylo.sor_nn8_null[,8,] + phylo.sor_nn8_null[,10,],
                        "alpha.min" = phylo.sor_nn8_null[,8,] + phylo.sor_nn8_null[,11,],
                        "diff.alpha" = phylo.sor_nn8_null[,10,] - phylo.sor_nn8_null[,11,],
                        "pnest" = phylo.sor_nn8_null[,6,]/phylo.sor_nn8_null[,7,],
                        along=2)
# get the mean and sd
phylo.sor_nn8_null_mean <- apply(phylo.sor_nn8_null, c(1,2), mean)
phylo.sor_nn8_null_sd <- apply(phylo.sor_nn8_null, c(1,2), sd)
phylo.sor_nn8_null_sd[,1:4] <- phylo.sor_nn8_null_mean[, 1:4]

## functional beta
# add maximum and minnimum alpha and their differences and percentage of nestedness
func.sor_nn8_null <- abind(func.sor_nn8_null, 
                        "alpha.max" = func.sor_nn8_null[,8,] + func.sor_nn8_null[,10,],
                        "alpha.min" = func.sor_nn8_null[,8,] + func.sor_nn8_null[,11,],
                        "diff.alpha" = func.sor_nn8_null[,10,] - func.sor_nn8_null[,11,],
                        "pnest" = func.sor_nn8_null[,6,]/func.sor_nn8_null[,7,],
                        along=2)
# get the mean and sd
func.sor_nn8_null_mean <- apply(func.sor_nn8_null, c(1, 2), mean)
func.sor_nn8_null_sd <- apply(func.sor_nn8_null, c(1, 2), sd)
func.sor_nn8_null_sd[,1:4] <- func.sor_nn8_null_mean[, 1:4]

save(phylo.sor_nn8_null_mean, phylo.sor_nn8_null_sd,
     func.sor_nn8_null_mean, func.sor_nn8_null_sd,
     file="intermediate_results/random_beta_phylo_func_nn8.RDATA")
