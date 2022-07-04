###################################################################################################
## Combine traits, phylogeny, remove species without distributions, generate dedrogram with traits####
###################################################################################################

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ape","picante", "fastcluster")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


# input distributions of all trees and traits and phylogeny
load("data/tree_pam/tree_pam6.RDATA")
tree <- read.tree("data/traits_phylogeny/TC_tree_v5.tre")
trait <- read_delim("data/traits_phylogeny/TC_species_8Traits_mean_final.csv", delim = ";")
spname <- read_delim("data/traits_phylogeny/TC_species.family.54020sp.csv", delim = ";")


# check species in traits and phylogeny can be matched
id <- tree$tip.label %in% trait[,1]
table(id)  #All is TRUE


##############
## prepare phylogenetic tree

# remove tips that did not included in species distributions
spp.keep <- tree$tip.label[tree$tip.label %in% tree_pam6$Spe]
tree <- keep.tip(tree, spp.keep)

# rootrize the tree
tree$root.edge <- 0
is.rooted(tree)

# convert the phylogeny into dichotomies to avoid errors in calculation of phylogenetic beta
tree <- multi2di(tree)

## check the phylogeny for angiosperm species
# angiosperm species
tree.tips <- tibble(species = tree$tip.label) %>%
  mutate(species2 = gsub("_"," ", species)) %>%
  left_join(spname, by = c("species2" = "species"))

table(tree.tips$class)

angio.spp <- tree.tips %>%
  filter(class == "Magnoliopsida") %>%
  pull(species)

# get the angiosperm phylogeny
tree.angio <- keep.tip(tree, angio.spp)

# plot the phylogeny for angiosperm species
pdf("data/traits_phylogeny/phylogeny_angiosperm.pdf", height=200, width=10)
plot(tree.angio, cex=0.05, edge.width=0.01)
dev.off()


###############
## prepare traits and dendrogram of traits

summary(trait)

# get the traits of species with distributions
trait <- trait %>%
  mutate(species2 = gsub("_"," ",species)) %>%
  left_join(spname, by = c("species2" = "species")) %>% 
   relocate(species, species2:class) %>%
   filter(species %in% tree_pam6$Species_name)

##Histogram of trait values
par(mfrow=c(3,3), mar=c(5,5,1,1))
hist(trait$leafN, xlab="Leaf N (mg/g)", main="", cex.lab=1.5, col="gray")
hist(trait$WDensity, xlab="Wood density (g/cm3)", main="", cex.lab=1.5, col="gray")
hist(trait$leafP, xlab="Leaf P (mg/g)", main="",cex.lab=1.5,col="gray")
hist(trait$LDMC, xlab="Leaf dry matter content (mg/g)", main="", cex.lab=1.5, col="gray")
hist(trait$VegHeight, xlab="Maximum height (m)", main="", cex.lab=1.5, col="gray")
hist(log10(trait$SeedDryMass), xlab="Log10(Seed dry mass (mg))", main="", cex.lab=1.5, col="gray")
hist(trait$SLA, xlab="Specific leaf area (mm2/mg)", main="", cex.lab=1.5, col="gray")
hist(log10(trait$LA), xlab="Log10(leaf area (mm2))", main="", cex.lab=1.5, col="gray")


## build trait dendrogram
# first, log-transform and then scale trait values
trait_scaled <- trait[, c(1, 7:14)]
trait_scaled[, c(2:9)] <- log(trait_scaled[, c(2:9)])
trait_scaled[, c(2:9)] <- scale(trait_scaled[, c(2:9)])

par(mfrow=c(3,3), mar=c(5,5,1,1))
for(i in 2:9){
  hist(trait_scaled[, i] %>% pull(), main = "")
}

# perform PCA for trait data
trait.pca <- princomp(trait_scaled[, 2:9])
summary(trait.pca)  #Cumulative Proportion  0.4272476 0.6773202 0.8363533 0.93461881 0.96439363
trait.pca$loading
trait.pca.score <- trait.pca$scores[, 1:4]
rownames(trait.pca.score) <- trait_scaled$species

# generate trait dendrogram, which will be used to calculate functional beta
trait.pca.dist <- dist(trait.pca.score)
trait.dendro <- fastcluster::hclust(trait.pca.dist)
trait.tree <- as.phylo(trait.dendro)

save(tree, trait, trait.tree, trait.pca.score, file = "data/traits_phylogeny/trait_phylogeny_final.RDATA")
