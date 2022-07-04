###################################################################################################
### apply the pineline that was used to prepare species distributions for ranges in alpha 6 
# to ranges that were estimated with alpha 2, 4, 10
###################################################################################################

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


##Tree distribution shape files
tree_range2 <- readOGR("data/tree_ranges/alpha2","TC_Sum2_raw")
tree_range4 <- readOGR("data/tree_ranges/alpha4","TC_Sum4_raw")
tree_range10 <- readOGR("data/tree_ranges/alpha10","TC_Sum10_raw")

# trait and species names
trait <- read_delim("data/traits_phylogeny/TC_species_8Traits_mean_final.csv", delim = ";")
spname <- read_delim("data/traits_phylogeny/TC_species.family.54020sp.csv", delim = ";")

# equal-area projection
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")

# global mainlands
land <- readOGR(dsn="data/global_maps/global_border_mainlands",layer="GSHHS_i_L1_simple")


##################
# For ranges producted with alpha = 2

# project tree ranges
summary(tree_range2)
colnames(tree_range2@data)<- "binomial"
tree_range2_bhm <- spTransform(tree_range2, behrmann)

# the world map and get the extent
data(wrld_simpl)
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world, behrmann)
ext <- extent(world)

# To keep running for spatial data with minor range exceedance of geographical coordinate range
print(bbox(tree_range2), digits=12)
set_ll_warn(FALSE) 
set_ll_TOL(0.2) 

# Convert species ranges into presence/absence matrix
tree_pam2 <- lets.presab(tree_range2_bhm, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], resol=200, crs = behrmann, crs.grid = behrmann)

# Set rownames of presence/absence matrix as the cell number of raster
rownames(tree_pam2[[1]]) <- cellFromXY(object=tree_pam2$Rich, xy=tree_pam2$Pre[,1:2])

save(tree_pam2, file="data/tree_pam/tree_pam2.RDATA")
# load("data/tree_pam/tree_pam2.RDATA")

# one cell in south America has very high richness, probably caused by the errors in distribution records
# the cell defined as outlier and its eight neighbour cells
id.out <- which(as.vector(tree_pam2$Rich > 3500))
row.nr <- rowFromCell(tree_pam2$Rich,id.out)
col.nr <- colFromCell(tree_pam2$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree_pam2$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

# remove species occurred in the cell but not in the eight neighbor cells
id.pre.out <- which(rownames(tree_pam2$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree_pam2$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree_pam2$Pre[id.pre.out.nb, -c(1:2)])
occ.out <- tree_pam2$Pre[id.pre.out, -c(1:2)]
id.spe.out <- which(occ.out.nb == 0 & occ.out == 1)
tree_pam2$Presence_and_Absence_Matrix[id.pre.out, id.spe.out+2] <- 0

#Update the richness of outlier cell
tree_pam2[[2]][id.out] <- sum(tree_pam2$Pre[id.pre.out,-c(1:2)])

# remove species without traits
id <- tree_pam2$Species_name %in% trait$species
table(id) #TRUE:46752; FALSE:15
spp <-  tree_pam2$Species_name[id]
tree_pam2 <- lets.subsetPAM(tree_pam2, names=spp, remove.cells = TRUE)


# Angiosperm species
angio_species <- tibble(species = tree_pam2$Spe) %>%
  mutate(species2 = gsub("_"," ", species)) %>%
  left_join(spname, by = c("species2" = "species")) %>%
  filter(class == "Magnoliopsida")

# only keep angiosperm species in the species presence/absence matrix
tree_pam2 <- lets.subsetPAM(tree_pam2, names = angio_species$species, remove.cells = TRUE)


## remove grid-cells with richness < 5
# the values of raster cells with richness <5 are set as NA 
tree_pam2$Richness_Raster[tree_pam2$Rich<5] <- NA 
table(!is.na(values(tree_pam2$Rich))) #TRUE:3098; FALSE:8734

## remove grid-cells that are not on mainland
land_bhm <- spTransform(land, behrmann)
tree_pam2$Richness_Raster <- mask(tree_pam2$Rich, land_bhm)  #based on mainland
table(!is.na(values(tree_pam2$Rich))) #FALSE:2319;TRUE:9513

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam2[[1]]) %in% which(!is.na(values(tree_pam2$Rich))) 
id_species <- colSums(tree_pam2[[1]][id_cell, -c(1:2)]) >0
tree_pam2$Presence_and_Absence_Matrix <- tree_pam2$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam2$Species_name <- colnames(tree_pam2[[1]])[-c(1:2)]

#Save data
save(tree_pam2, file = "data/tree_pam/tree_pam2_final.RDATA")


##################
# For ranges producted with alpha = 4

# project tree ranges
summary(tree_range4)
colnames(tree_range4@data)<- "binomial"
tree_range4_bhm <- spTransform(tree_range4, behrmann)

# the world map and get the extent
data(wrld_simpl)
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world, behrmann)
ext <- extent(world)

# To keep running for spatial data with minor range exceedance of geographical coordinate range
print(bbox(tree_range4), digits=12)
set_ll_warn(FALSE) 
set_ll_TOL(0.2) 

# Convert species ranges into presence/absence matrix
tree_pam4 <- lets.presab(tree_range4_bhm, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], resol=200, crs = behrmann, crs.grid = behrmann)

# Set rownames of presence/absence matrix as the cell number of raster
rownames(tree_pam4[[1]]) <- cellFromXY(object=tree_pam4$Rich, xy=tree_pam4$Pre[,1:2])

save(tree_pam4, file="data/tree_pam/tree_pam4.RDATA")
# load("data/tree_pam/tree_pam4.RDATA")

# one cell in south America has very high richness, probably caused by the errors in distribution records
# the cell defined as outlier and its eight neighbour cells
id.out <- which(as.vector(tree_pam4$Rich > 3500))
row.nr <- rowFromCell(tree_pam4$Rich,id.out)
col.nr <- colFromCell(tree_pam4$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree_pam4$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

# remove species occurred in the cell but not in the eight neighbor cells
id.pre.out <- which(rownames(tree_pam4$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree_pam4$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree_pam4$Pre[id.pre.out.nb, -c(1:2)])
occ.out <- tree_pam4$Pre[id.pre.out, -c(1:2)]
id.spe.out <- which(occ.out.nb == 0 & occ.out == 1)
tree_pam4$Presence_and_Absence_Matrix[id.pre.out, id.spe.out+2] <- 0

#Update the richness of outlier cell
tree_pam4[[2]][id.out] <- sum(tree_pam4$Pre[id.pre.out,-c(1:2)])

# remove species without traits
id <- tree_pam4$Species_name %in% trait$species
table(id) #TRUE:46752; FALSE:15
spp <-  tree_pam4$Species_name[id]
tree_pam4 <- lets.subsetPAM(tree_pam4, names=spp, remove.cells = TRUE)


# Angiosperm species
angio_species <- tibble(species = tree_pam4$Spe) %>%
  mutate(species2 = gsub("_"," ", species)) %>%
  left_join(spname, by = c("species2" = "species")) %>%
  filter(class == "Magnoliopsida")

# only keep angiosperm species in the species presence/absence matrix
tree_pam4 <- lets.subsetPAM(tree_pam4, names = angio_species$species, remove.cells = TRUE)


## remove grid-cells with richness < 5
# the values of raster cells with richness <5 are set as NA 
tree_pam4$Richness_Raster[tree_pam4$Rich<5] <- NA 
table(!is.na(values(tree_pam4$Rich))) #TRUE:3098; FALSE:8734

## remove grid-cells that are not on mainland
land_bhm <- spTransform(land, behrmann)
tree_pam4$Richness_Raster <- mask(tree_pam4$Rich, land_bhm)  #based on mainland
table(!is.na(values(tree_pam4$Rich))) #FALSE:2319;TRUE:9513

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam4[[1]]) %in% which(!is.na(values(tree_pam4$Rich))) 
id_species <- colSums(tree_pam4[[1]][id_cell, -c(1:2)]) >0
tree_pam4$Presence_and_Absence_Matrix <- tree_pam4$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam4$Species_name <- colnames(tree_pam4[[1]])[-c(1:2)]

#Save data
save(tree_pam4, file = "data/tree_pam/tree_pam4_final.RDATA")


##################
# For ranges producted with alpha = 10

# project tree ranges
summary(tree_range10)
colnames(tree_range10@data)<- "binomial"
tree_range10_bhm <- spTransform(tree_range10, behrmann)

# the world map and get the extent
data(wrld_simpl)
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world, behrmann)
ext <- extent(world)

# To keep running for spatial data with minor range exceedance of geographical coordinate range
print(bbox(tree_range10), digits=12)
set_ll_warn(FALSE) 
set_ll_TOL(0.2) 

# Convert species ranges into presence/absence matrix
tree_pam10 <- lets.presab(tree_range10_bhm, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], resol=200, crs = behrmann, crs.grid = behrmann)

# Set rownames of presence/absence matrix as the cell number of raster
rownames(tree_pam10[[1]]) <- cellFromXY(object=tree_pam10$Rich, xy=tree_pam10$Pre[,1:2])

save(tree_pam10, file="data/tree_pam/tree_pam10.RDATA")
# load("data/tree_pam/tree_pam10.RDATA")

# one cell in south America has very high richness, probably caused by the errors in distribution records
# the cell defined as outlier and its eight neighbour cells
id.out <- which(as.vector(tree_pam10$Rich > 3500))
row.nr <- rowFromCell(tree_pam10$Rich,id.out)
col.nr <- colFromCell(tree_pam10$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree_pam10$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

# remove species occurred in the cell but not in the eight neighbor cells
id.pre.out <- which(rownames(tree_pam10$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree_pam10$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree_pam10$Pre[id.pre.out.nb, -c(1:2)])
occ.out <- tree_pam10$Pre[id.pre.out, -c(1:2)]
id.spe.out <- which(occ.out.nb == 0 & occ.out == 1)
tree_pam10$Presence_and_Absence_Matrix[id.pre.out, id.spe.out+2] <- 0

#Update the richness of outlier cell
tree_pam10[[2]][id.out] <- sum(tree_pam10$Pre[id.pre.out,-c(1:2)])

# remove species without traits
id <- tree_pam10$Species_name %in% trait$species
table(id) #TRUE:46752; FALSE:15
spp <-  tree_pam10$Species_name[id]
tree_pam10 <- lets.subsetPAM(tree_pam10, names=spp, remove.cells = TRUE)


# Angiosperm species
angio_species <- tibble(species = tree_pam10$Spe) %>%
  mutate(species2 = gsub("_"," ", species)) %>%
  left_join(spname, by = c("species2" = "species")) %>%
  filter(class == "Magnoliopsida")

# only keep angiosperm species in the species presence/absence matrix
tree_pam10 <- lets.subsetPAM(tree_pam10, names = angio_species$species, remove.cells = TRUE)


## remove grid-cells with richness < 5
# the values of raster cells with richness <5 are set as NA 
tree_pam10$Richness_Raster[tree_pam10$Rich<5] <- NA 
table(!is.na(values(tree_pam10$Rich))) #TRUE:3098; FALSE:8734

## remove grid-cells that are not on mainland
land_bhm <- spTransform(land, behrmann)
tree_pam10$Richness_Raster <- mask(tree_pam10$Rich, land_bhm)  #based on mainland
table(!is.na(values(tree_pam10$Rich))) #FALSE:2319;TRUE:9513

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam10[[1]]) %in% which(!is.na(values(tree_pam10$Rich))) 
id_species <- colSums(tree_pam10[[1]][id_cell, -c(1:2)]) >0
tree_pam10$Presence_and_Absence_Matrix <- tree_pam10$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam10$Species_name <- colnames(tree_pam10[[1]])[-c(1:2)]

#Save data
save(tree_pam10, file = "data/tree_pam/tree_pam10_final.RDATA")

