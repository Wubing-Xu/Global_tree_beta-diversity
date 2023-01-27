###################################################################################################
## Calculate the number of species, genus and family that have observational trait data and  their completeness 
## calculate trait coverage (proportion of species with at least one to eight traits within each grid cell)

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse", "rgeos", "rgdal", "raster", "maptools", "sf", "letsR")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


# input observational and imputed trait data and taxonomic names
trait_mean <- read_delim("data/traits_phylogeny/traitsOutMatrix_Observations2Traits_v5.mean.csv", delim = ";")[, -1]
trait_imputed <- read_delim("data/traits_phylogeny/TC_species_8Traits_mean_final.csv", delim = ";")
spname <- read_delim("data/traits_phylogeny/TC_species.family.54020sp.csv", delim = ";")

# add genus, family, Order, class for species names, and tranform the format
trait_mean <- trait_mean %>%
  left_join(spname, by = "species") %>% 
  relocate(species, genus:class) %>%
  pivot_longer(cols = leafN:WoodDM, names_to = "trait_name", values_to = "trait_mean") %>%
  dplyr::filter(!is.na(trait_mean))

trait_imputed <- trait_imputed %>%
  mutate(species = gsub("_", " ", species)) %>%
  left_join(spname, by = "species") %>% 
  relocate(species, genus:class) %>%
  pivot_longer(cols = leafN:LA, names_to = "trait_name", values_to = "trait_values") 

# get the trait values for species with observational data from the two tables (mean and imputed)
trait_obs <- trait_mean %>% 
  dplyr::filter(trait_name %in% unique(trait_imputed$trait_name)) %>%
  left_join(trait_imputed)


# compare trait values from two tables
# they were slightly different (particular LDMC, VegHeight); 
# this is because the observed was the mean value, while the one in imputed table consider trait distributions (more reliable; used for furhter)
ggplot(trait_obs %>% filter(trait_name == "LA")) +
  geom_point(aes(x = trait_mean, y = trait_values)) +
  geom_abline() +
  scale_x_log10() +
  scale_y_log10()

####################
# calculate the number of species, genus and family that have observational trait data and their completeness 

# total number of tree species, genus and family
total_taxa <- spname %>% 
  summarise(nspecies = n_distinct(species),
            ngenus = n_distinct(genus),
            nfamily = n_distinct(family))

# total number of species, genus and family that have observational trait data
trait_total_taxa <- trait_mean %>% 
  summarise(nspecies = n_distinct(species[!is.na(trait_mean)]),
            ngenus = n_distinct(genus[!is.na(trait_mean)]),
            nfamily = n_distinct(family[!is.na(trait_mean)])) %>%
  mutate(trait_name = "all")

# calculate the number of species, genus and family that have observational trait data for each trait and 
# their completeness (proportion of species of total taxa)
trait_completeness <- trait_mean %>% 
  group_by(trait_name) %>%
  summarise(nspecies = n_distinct(species[!is.na(trait_mean)]),
            ngenus = n_distinct(genus[!is.na(trait_mean)]),
            nfamily = n_distinct(family[!is.na(trait_mean)])) %>% 
  bind_rows(trait_total_taxa) %>%
  mutate(complete_species = round(100*nspecies/total_taxa$nspecies, 1),
         complete_genus = round(100*ngenus/total_taxa$ngenus, 1),
         complete_family = round(100*nfamily/total_taxa$nfamily, 1))

write_csv(trait_completeness, file = "results/trait_completeness.csv")



####################
## calculate trait coverage (proportion of species with at least one to eight traits within each grid cell)

load("data/tree_pam/tree_pam6_final.RDATA")
trait <- read_delim("data/traits_phylogeny/TC_species_8Traits_mean_final.csv", delim = ";")

# change the species names format to match other files 
trait_obs <- trait_obs %>%
  rename(species2 = species) %>%
  mutate(species = gsub(" ","_",species2)) %>%
  relocate(species, species2)

# the number of observational traits for each species
spname_ntraits <- trait_obs %>%
  group_by(species) %>%
  summarise(ntraits_obs = n_distinct(trait_name))

trait <- trait %>% 
  left_join(spname_ntraits) %>%
   mutate(ntraits_obs = ifelse(is.na(ntraits_obs), 0, ntraits_obs))

# keep only species with distributions: a total of 43635 species
trait <- trait %>%
  filter(species %in% tree_pam6$Species_name)

# the number of species with at least one, three or all traits
trait %>% filter(ntraits_obs >= 1) #10,635 species
trait %>% filter(ntraits_obs >= 2) # 4999 species
trait %>% filter(ntraits_obs >= 3) # 3230 species
trait %>% filter(ntraits_obs >= 4) # 2177 species
trait %>% filter(ntraits_obs >= 5) # 1264 species
trait %>% filter(ntraits_obs >= 6) # 3230 species
trait %>% filter(ntraits_obs >= 7) # 646 species
trait %>% filter(ntraits_obs == 8) # only 114 species


## update distribution data: keep only species with at least one to eight observational traits
species_1trait <- trait %>% filter(ntraits_obs >= 1) %>% pull(species)
species_2trait <- trait %>% filter(ntraits_obs >= 2) %>% pull(species)
species_3trait <- trait %>% filter(ntraits_obs >= 3) %>% pull(species)
species_4trait <- trait %>% filter(ntraits_obs >= 4) %>% pull(species)
species_5trait <- trait %>% filter(ntraits_obs >= 5) %>% pull(species)
species_6trait <- trait %>% filter(ntraits_obs >= 6) %>% pull(species)
species_7trait <- trait %>% filter(ntraits_obs >= 7) %>% pull(species)
species_8trait <- trait %>% filter(ntraits_obs >= 8) %>% pull(species)

tree_pam1trait <- lets.subsetPAM(tree_pam6, names = species_1trait, remove.cells = TRUE)
tree_pam2trait <- lets.subsetPAM(tree_pam6, names = species_2trait, remove.cells = TRUE)
tree_pam3trait <- lets.subsetPAM(tree_pam6, names = species_3trait, remove.cells = TRUE)
tree_pam4trait <- lets.subsetPAM(tree_pam6, names = species_4trait, remove.cells = TRUE)
tree_pam5trait <- lets.subsetPAM(tree_pam6, names = species_5trait, remove.cells = TRUE)
tree_pam6trait <- lets.subsetPAM(tree_pam6, names = species_6trait, remove.cells = TRUE)
tree_pam7trait <- lets.subsetPAM(tree_pam6, names = species_7trait, remove.cells = TRUE)
tree_pam8trait <- lets.subsetPAM(tree_pam6, names = species_8trait, remove.cells = TRUE)

# calculate trait coverage: number and proportion of species with trait data
tree_sprich_all <- tibble(ID = as.numeric(rownames(tree_pam6[[1]])), sprich_all = apply(tree_pam6[[1]][, -c(1:2)], 1, sum))
tree_sprich_1trait <- tibble(ID = as.numeric(rownames(tree_pam1trait[[1]])), sprich_1trait = apply(tree_pam1trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_2trait <- tibble(ID = as.numeric(rownames(tree_pam2trait[[1]])), sprich_2trait = apply(tree_pam2trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_3trait <- tibble(ID = as.numeric(rownames(tree_pam3trait[[1]])), sprich_3trait = apply(tree_pam3trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_4trait <- tibble(ID = as.numeric(rownames(tree_pam4trait[[1]])), sprich_4trait = apply(tree_pam4trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_5trait <- tibble(ID = as.numeric(rownames(tree_pam5trait[[1]])), sprich_5trait = apply(tree_pam5trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_6trait <- tibble(ID = as.numeric(rownames(tree_pam6trait[[1]])), sprich_6trait = apply(tree_pam6trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_7trait <- tibble(ID = as.numeric(rownames(tree_pam7trait[[1]])), sprich_7trait = apply(tree_pam7trait[[1]][, -c(1:2)], 1, sum))
tree_sprich_8trait <- tibble(ID = as.numeric(rownames(tree_pam8trait[[1]])), sprich_8trait = apply(tree_pam8trait[[1]][, -c(1:2)], 1, sum))

trait_coverage <- tree_sprich_all %>%
  left_join(tree_sprich_1trait) %>%
  left_join(tree_sprich_2trait) %>%
  left_join(tree_sprich_3trait) %>%
  left_join(tree_sprich_4trait) %>%
  left_join(tree_sprich_5trait) %>%
  left_join(tree_sprich_6trait) %>%
  left_join(tree_sprich_7trait) %>%
  left_join(tree_sprich_8trait) %>%
  mutate(sprich_1trait = ifelse(is.na(sprich_1trait), 0, sprich_1trait),
         sprich_2trait = ifelse(is.na(sprich_2trait), 0, sprich_2trait),
         sprich_3trait = ifelse(is.na(sprich_3trait), 0, sprich_3trait),
         sprich_4trait = ifelse(is.na(sprich_4trait), 0, sprich_4trait),
         sprich_5trait = ifelse(is.na(sprich_5trait), 0, sprich_5trait),
         sprich_6trait = ifelse(is.na(sprich_6trait), 0, sprich_6trait),
         sprich_7trait = ifelse(is.na(sprich_7trait), 0, sprich_7trait),
         sprich_8trait = ifelse(is.na(sprich_8trait), 0, sprich_8trait),
         prop_spp_1trait = sprich_1trait/sprich_all,
         prop_spp_2trait = sprich_2trait/sprich_all,
         prop_spp_3trait = sprich_3trait/sprich_all,
         prop_spp_4trait = sprich_4trait/sprich_all,
         prop_spp_5trait = sprich_5trait/sprich_all,
         prop_spp_6trait = sprich_6trait/sprich_all,
         prop_spp_7trait = sprich_7trait/sprich_all,
         prop_spp_8trait = sprich_8trait/sprich_all) %>%
  dplyr::select(ID, sprich_all, prop_spp_1trait:prop_spp_8trait)

summary(trait_coverage$prop_spp_1trait)
summary(trait_coverage$prop_spp_3trait)
summary(trait_coverage$prop_spp_5trait)
summary(trait_coverage$prop_spp_7trait)


########
##update tree_pam files for calculating beta-diversity

## remove grid-cells with richness < 5
table(values(tree_pam1trait$Rich) >= 5)
table(values(tree_pam3trait$Rich) >= 5)
table(values(tree_pam5trait$Rich) >= 5)

#### for tree_pam1trait
# the values of raster cells with richness <5 are set as NA 
tree_pam1trait$Richness_Raster[tree_pam1trait$Rich<5] <- NA 

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam1trait[[1]]) %in% which(!is.na(values(tree_pam1trait$Rich))) 
id_species <- colSums(tree_pam1trait[[1]][id_cell, -c(1:2)]) >0
tree_pam1trait$Presence_and_Absence_Matrix <- tree_pam1trait$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam1trait$Species_name <- colnames(tree_pam1trait[[1]])[-c(1:2)]


#### for tree_pam3trait
# the values of raster cells with richness <5 are set as NA 
tree_pam3trait$Richness_Raster[tree_pam3trait$Rich<5] <- NA 

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam3trait[[1]]) %in% which(!is.na(values(tree_pam3trait$Rich))) 
id_species <- colSums(tree_pam3trait[[1]][id_cell, -c(1:2)]) >0
tree_pam3trait$Presence_and_Absence_Matrix <- tree_pam3trait$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam3trait$Species_name <- colnames(tree_pam3trait[[1]])[-c(1:2)]


#### for tree_pam5trait
# the values of raster cells with richness <5 are set as NA 
tree_pam5trait$Richness_Raster[tree_pam5trait$Rich<5] <- NA 

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam5trait[[1]]) %in% which(!is.na(values(tree_pam5trait$Rich))) 
id_species <- colSums(tree_pam5trait[[1]][id_cell, -c(1:2)]) >0
tree_pam5trait$Presence_and_Absence_Matrix <- tree_pam5trait$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam5trait$Species_name <- colnames(tree_pam5trait[[1]])[-c(1:2)]

# save files
save(tree_pam1trait, tree_pam3trait, tree_pam5trait, 
     file = "data/traits_sensitivity/tree_pam6_SpeciesWithTraits.RDATA")

save(trait_obs, trait, trait_coverage, 
     file = "data/traits_sensitivity/traits_coverage_sensitivity.RDATA")

