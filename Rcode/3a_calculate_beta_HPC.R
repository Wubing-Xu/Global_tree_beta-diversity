############################################################################################################################
## Calculate taxonomic/phylogenetic/functional beta and its turnvoer and nesedness and proportion of nestedness 
## among each focal cell and its 13 (or 5) of 24 (or 9) nearest neighboring cells using multiple-site sorensen dissimilarity
## Perform these calculations using different range estimates and size of moving windows 
############################################################################################################################

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

# load the self_defined function for calculation of beta-diversity
source("Rcode/functions/get_beta.R")

load("data/tree_pam/tree_pam6_final.RDATA")
load("data/tree_pam/tree_pam2_final.RDATA")
load("data/tree_pam/tree_pam4_final.RDATA")
load("data/tree_pam/tree_pam10_final.RDATA")
load("data/traits_phylogeny/trait_phylogeny_final.RDATA")


# set parallel computing
t1 <- Sys.time()
print("n_cores")
print(detectCores()) #the number of cores
no_cores <- 35   
cl <- makeCluster(no_cores, type="PSOCK") # Initiate cluster
clusterExport(cl, varlist = c("tree_pam6", "tree_pam2", "tree_pam4", "tree_pam10", "tree", "trait.tree", "get_beta"))
clusterEvalQ(cl, c(library(betapart,), library(spdep, lib.loc = "/gpfs0/home/wubing/R/library"), 
			library(picante), library(ape), library(phytools), library(abind)))

## using species distributions estimated with alpha hulls in 6
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
spe.sor <- parLapply(cl,1:200, get_beta_tempo)
spe.sor <- sapply(spe.sor, as.matrix, simplify="array")
spe.sor <- abind(spe.sor, "pnest" = spe.sor[,6,]/spe.sor[,7,], along=2)
spe.sor <- apply(spe.sor, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
phylo.sor <- parLapply(cl,1:200, get_beta_tempo)
phylo.sor <- sapply(phylo.sor, as.matrix, simplify="array")
phylo.sor <- abind(phylo.sor, "pnest" = phylo.sor[,6,]/phylo.sor[,7,], along=2)
phylo.sor <- apply(phylo.sor, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
func.sor <- parLapply(cl,1:200, get_beta_tempo)
func.sor <- sapply(func.sor, as.matrix, simplify="array")
func.sor <- abind(func.sor, "pnest" = func.sor[,6,]/func.sor[,7,], along=2)
func.sor <- apply(func.sor, c(1, 2), mean)

save(spe.sor, phylo.sor, func.sor, 
     file="intermediate_results/observed_beta_sorensen_alpha6.RDATA")
t2 <- Sys.time()
t2 - t1


# calculate beta-diversity among focal cell and its 5 of 8 nearest neighboring cells (3*3 cells)
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=5, d2=300, ncell_resample=5), silent=TRUE)
spe.sor_nn8 <- parLapply(cl,1:200, get_beta_tempo)
spe.sor_nn8 <- sapply(spe.sor_nn8, as.matrix, simplify="array")
spe.sor_nn8 <- abind(spe.sor_nn8, "pnest" = spe.sor_nn8[,6,]/spe.sor_nn8[,7,], along=2)
spe.sor_nn8 <- apply(spe.sor_nn8, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=5, d2=300, ncell_resample=5), silent=TRUE)
phylo.sor_nn8 <- parLapply(cl,1:200, get_beta_tempo)
phylo.sor_nn8 <- sapply(phylo.sor_nn8, as.matrix, simplify="array")
phylo.sor_nn8 <- abind(phylo.sor_nn8, "pnest" = phylo.sor_nn8[,6,]/phylo.sor_nn8[,7,], along=2)
phylo.sor_nn8 <- apply(phylo.sor_nn8, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam6$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=5, d2=300, ncell_resample=5), silent=TRUE)
func.sor_nn8 <- parLapply(cl,1:200, get_beta_tempo)
func.sor_nn8 <- sapply(func.sor_nn8, as.matrix, simplify="array")
func.sor_nn8 <- abind(func.sor_nn8, "pnest" = func.sor_nn8[,6,]/func.sor_nn8[,7,], along=2)
func.sor_nn8 <- apply(func.sor_nn8, c(1, 2), mean)

save(spe.sor_nn8, phylo.sor_nn8, func.sor_nn8, 
     file="intermediate_results/observed_beta_sorensen_nn8.RDATA")
t3 <- Sys.time()
t3 - t2


# using species distributions estimated with alpha hulls in 2
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam2$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
spe.sor_pam2 <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam2 <- sapply(spe.sor_pam2, as.matrix, simplify="array")
spe.sor_pam2 <- abind(spe.sor_pam2, "pnest" = spe.sor_pam2[,6,]/spe.sor_pam2[,7,], along=2)
spe.sor_pam2 <- apply(spe.sor_pam2, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam2$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
phylo.sor_pam2 <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam2 <- sapply(phylo.sor_pam2, as.matrix, simplify="array")
phylo.sor_pam2 <- abind(phylo.sor_pam2, "pnest" = phylo.sor_pam2[,6,]/phylo.sor_pam2[,7,], along=2)
phylo.sor_pam2 <- apply(phylo.sor_pam2, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam2$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
func.sor_pam2 <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam2 <- sapply(func.sor_pam2, as.matrix, simplify="array")
func.sor_pam2 <- abind(func.sor_pam2, "pnest" = func.sor_pam2[,6,]/func.sor_pam2[,7,], along=2)
func.sor_pam2 <- apply(func.sor_pam2, c(1, 2), mean)

save(spe.sor_pam2, phylo.sor_pam2, func.sor_pam2, 
     file="intermediate_results/observed_beta_sorensen_alpha2.RDATA")
t4 <- Sys.time()
t4 - t3


# using species distributions estimated with alpha hulls in 4
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam4$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
spe.sor_pam4 <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam4 <- sapply(spe.sor_pam4, as.matrix, simplify="array")
spe.sor_pam4 <- abind(spe.sor_pam4, "pnest" = spe.sor_pam4[,6,]/spe.sor_pam4[,7,], along=2)
spe.sor_pam4 <- apply(spe.sor_pam4, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam4$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
phylo.sor_pam4 <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam4 <- sapply(phylo.sor_pam4, as.matrix, simplify="array")
phylo.sor_pam4 <- abind(phylo.sor_pam4, "pnest" = phylo.sor_pam4[,6,]/phylo.sor_pam4[,7,], along=2)
phylo.sor_pam4 <- apply(phylo.sor_pam4, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam4$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
func.sor_pam4 <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam4 <- sapply(func.sor_pam4, as.matrix, simplify="array")
func.sor_pam4 <- abind(func.sor_pam4, "pnest" = func.sor_pam4[,6,]/func.sor_pam4[,7,], along=2)
func.sor_pam4 <- apply(func.sor_pam4, c(1, 2), mean)

save(spe.sor_pam4, phylo.sor_pam4, func.sor_pam4, 
     file="intermediate_results/observed_beta_sorensen_alpha4.RDATA")
t5 <- Sys.time()
t5 - t4


# using species distributions estimated with alpha hulls in 10
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam10$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
spe.sor_pam10 <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam10 <- sapply(spe.sor_pam10, as.matrix, simplify="array")
spe.sor_pam10 <- abind(spe.sor_pam10, "pnest" = spe.sor_pam10[,6,]/spe.sor_pam10[,7,], along=2)
spe.sor_pam10 <- apply(spe.sor_pam10, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam10$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
phylo.sor_pam10 <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam10 <- sapply(phylo.sor_pam10, as.matrix, simplify="array")
phylo.sor_pam10 <- abind(phylo.sor_pam10, "pnest" = phylo.sor_pam10[,6,]/phylo.sor_pam10[,7,], along=2)
phylo.sor_pam10 <- apply(phylo.sor_pam10, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam10$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample=13), silent=TRUE)
func.sor_pam10 <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam10 <- sapply(func.sor_pam10, as.matrix, simplify="array")
func.sor_pam10 <- abind(func.sor_pam10, "pnest" = func.sor_pam10[,6,]/func.sor_pam10[,7,], along=2)
func.sor_pam10 <- apply(func.sor_pam10, c(1, 2), mean)

save(spe.sor_pam10, phylo.sor_pam10, func.sor_pam10, 
     file="intermediate_results/observed_beta_sorensen_alpha10.RDATA")
t6 <- Sys.time()
t6 - t5

save(spe.sor, phylo.sor, func.sor, 
     spe.sor_pam2, phylo.sor_pam2, func.sor_pam2,
     spe.sor_pam4, phylo.sor_pam4, func.sor_pam4,
     spe.sor_pam10, phylo.sor_pam10, func.sor_pam10,
     spe.sor_nn8, phylo.sor_nn8, func.sor_nn8,
     file="intermediate_results/observed_beta_sorensen.RDATA")

stopCluster(cl)
