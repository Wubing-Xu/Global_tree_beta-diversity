############################################################################################################################
## Calculate taxonomic/phylogenetic/functional beta and its turnvoer and nesedness and proportion of nestedness using species with observational traits
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

load("data/traits_sensitivity/tree_pam6_SpeciesWithTraits.RDATA")
load("data/traits_phylogeny/trait_phylogeny_final.RDATA")


# set parallel computing
t3 <- Sys.time()
print("n_cores")
print(detectCores()) #the number of cores
no_cores <- 25  
cl <- makeCluster(no_cores, type="PSOCK") # Initiate cluster
clusterExport(cl, varlist = c("tree_pam1trait", "tree_pam3trait", "tree_pam5trait", "tree", "trait.tree", "get_beta"))
clusterEvalQ(cl, c(library(betapart,), library(spdep, lib.loc = "/gpfs0/home/wubing/R/library"), 
			library(picante), library(ape), library(phytools), library(abind)))

# using only species with at least one measured trait
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam1trait$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
spe.sor_pam1trait <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam1trait <- sapply(spe.sor_pam1trait, as.matrix, simplify="array")
spe.sor_pam1trait <- abind(spe.sor_pam1trait, "pnest" = spe.sor_pam1trait[,6,]/spe.sor_pam1trait[,7,], along=2)
spe.sor_pam1trait <- apply(spe.sor_pam1trait, c(1, 2), mean)

t3 

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam1trait$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
phylo.sor_pam1trait <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam1trait <- sapply(phylo.sor_pam1trait, as.matrix, simplify="array")
phylo.sor_pam1trait <- abind(phylo.sor_pam1trait, "pnest" = phylo.sor_pam1trait[,6,]/phylo.sor_pam1trait[,7,], along=2)
phylo.sor_pam1trait <- apply(phylo.sor_pam1trait, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam1trait$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
func.sor_pam1trait <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam1trait <- sapply(func.sor_pam1trait, as.matrix, simplify="array")
func.sor_pam1trait <- abind(func.sor_pam1trait, "pnest" = func.sor_pam1trait[,6,]/func.sor_pam1trait[,7,], along=2)
func.sor_pam1trait <- apply(func.sor_pam1trait, c(1, 2), mean)

save(spe.sor_pam1trait, phylo.sor_pam1trait, func.sor_pam1trait, 
     file="intermediate_results/observed_beta_sorensen_SpeciesWith1trait.RDATA")
t4 <- Sys.time()
t4 - t3


# using only species with at least three measured traits
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam3trait$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
spe.sor_pam3trait <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam3trait <- sapply(spe.sor_pam3trait, as.matrix, simplify="array")
spe.sor_pam3trait <- abind(spe.sor_pam3trait, "pnest" = spe.sor_pam3trait[,6,]/spe.sor_pam3trait[,7,], along=2)
spe.sor_pam3trait <- apply(spe.sor_pam3trait, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam3trait$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
phylo.sor_pam3trait <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam3trait <- sapply(phylo.sor_pam3trait, as.matrix, simplify="array")
phylo.sor_pam3trait <- abind(phylo.sor_pam3trait, "pnest" = phylo.sor_pam3trait[,6,]/phylo.sor_pam3trait[,7,], along=2)
phylo.sor_pam3trait <- apply(phylo.sor_pam3trait, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam3trait$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
func.sor_pam3trait <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam3trait <- sapply(func.sor_pam3trait, as.matrix, simplify="array")
func.sor_pam3trait <- abind(func.sor_pam3trait, "pnest" = func.sor_pam3trait[,6,]/func.sor_pam3trait[,7,], along=2)
func.sor_pam3trait <- apply(func.sor_pam3trait, c(1, 2), mean)

save(spe.sor_pam3trait, phylo.sor_pam3trait, func.sor_pam3trait, 
     file="intermediate_results/observed_beta_sorensen_SpeciesWith3trait.RDATA")
t5 <- Sys.time()
t5 - t4


# using only species with at least five measured traits
# taxonomic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam5trait$Pre, method="sorensen", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
spe.sor_pam5trait <- parLapply(cl,1:100, get_beta_tempo)
spe.sor_pam5trait <- sapply(spe.sor_pam5trait, as.matrix, simplify="array")
spe.sor_pam5trait <- abind(spe.sor_pam5trait, "pnest" = spe.sor_pam5trait[,6,]/spe.sor_pam5trait[,7,], along=2)
spe.sor_pam5trait <- apply(spe.sor_pam5trait, c(1, 2), mean)

# phylogenetic beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam5trait$Pre,tree=tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
phylo.sor_pam5trait <- parLapply(cl,1:100, get_beta_tempo)
phylo.sor_pam5trait <- sapply(phylo.sor_pam5trait, as.matrix, simplify="array")
phylo.sor_pam5trait <- abind(phylo.sor_pam5trait, "pnest" = phylo.sor_pam5trait[,6,]/phylo.sor_pam5trait[,7,], along=2)
phylo.sor_pam5trait <- apply(phylo.sor_pam5trait, c(1, 2), mean)

# functional beta
get_beta_tempo <- function(i) try(get_beta(comdata=tree_pam5trait$Pre,tree=trait.tree, method="phylosor", beta.multi=TRUE, alpha=TRUE, n=13, d2=570, ncell_resample = 13), silent=TRUE)
func.sor_pam5trait <- parLapply(cl,1:100, get_beta_tempo)
func.sor_pam5trait <- sapply(func.sor_pam5trait, as.matrix, simplify="array")
func.sor_pam5trait <- abind(func.sor_pam5trait, "pnest" = func.sor_pam5trait[,6,]/func.sor_pam5trait[,7,], along=2)
func.sor_pam5trait <- apply(func.sor_pam5trait, c(1, 2), mean)

save(spe.sor_pam5trait, phylo.sor_pam5trait, func.sor_pam5trait, 
     file="intermediate_results/observed_beta_sorensen_SpeciesWith5trait.RDATA")
t6 <- Sys.time()
t6 - t5


save(spe.sor_pam1trait, phylo.sor_pam1trait, func.sor_pam1trait, 
     spe.sor_pam3trait, phylo.sor_pam3trait, func.sor_pam3trait, 
     spe.sor_pam5trait, phylo.sor_pam5trait, func.sor_pam5trait,
     file = "intermediate_results/observed_beta_sorensen_SpeciesWithTraits.RDATA")

stopCluster(cl)
