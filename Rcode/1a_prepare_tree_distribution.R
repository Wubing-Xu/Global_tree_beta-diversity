###################################################################################################
###Transform tree ranges to 'PAM' format with site-species matrix, 
# update data in one outlier cell in South America, keep species with trait data and are angiosperm; 
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
tree_range6 <- readOGR("data/tree_ranges/alpha6","TC_Sum6_raw")

# trait and species names
trait <- read_delim("data/traits_phylogeny/TC_species_8Traits_mean_final.csv", delim = ";")
spname <- read_delim("data/traits_phylogeny/TC_species.family.54020sp.csv", delim = ";")

# equal-area projection
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")

# global mainlands
land <- readOGR(dsn="data/global_maps/global_border_mainlands",layer="GSHHS_i_L1_simple")


# project tree ranges
summary(tree_range6)
colnames(tree_range6@data)<- "binomial"
tree_range6_bhm <- spTransform(tree_range6, behrmann)

# the world map and get the extent
data(wrld_simpl)
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world, behrmann)
ext <- extent(world)

# To keep running for spatial data with minor range exceedance of geographical coordinate range
print(bbox(tree_range6), digits=12)
set_ll_warn(FALSE) 
set_ll_TOL(0.2) 

# Convert species ranges into presence/absence matrix
tree_pam6 <- lets.presab(tree_range6_bhm, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], resol=200, crs = behrmann, crs.grid = behrmann)

# Set rownames of presence/absence matrix as the cell number of raster
rownames(tree_pam6[[1]]) <- cellFromXY(object=tree_pam6$Rich, xy=tree_pam6$Pre[,1:2])

save(tree_pam6, file="data/tree_pam/tree_pam6.RDATA")
# load("data/tree_pam/tree_pam6.RDATA")


# one cell in south America has very high richness, probably caused by the errors in distribution records
# the cell defined as outlier and its eight neighbour cells
id.out <- which(as.vector(tree_pam6$Rich > 3500))
row.nr <- rowFromCell(tree_pam6$Rich,id.out)
col.nr <- colFromCell(tree_pam6$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree_pam6$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

# remove species occurred in the cell but not in the eight neighbor cells
id.pre.out <- which(rownames(tree_pam6$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree_pam6$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree_pam6$Pre[id.pre.out.nb, -c(1:2)])
occ.out <- tree_pam6$Pre[id.pre.out, -c(1:2)]
id.spe.out <- which(occ.out.nb == 0 & occ.out == 1)
tree_pam6$Presence_and_Absence_Matrix[id.pre.out, id.spe.out+2] <- 0

#Update the richness of outlier cell
tree_pam6[[2]][id.out] <- sum(tree_pam6$Pre[id.pre.out,-c(1:2)])

# remove species without traits
id <- tree_pam6$Species_name %in% trait$species
table(id) #TRUE:46752; FALSE:15
spp <-  tree_pam6$Species_name[id]
tree_pam6 <- lets.subsetPAM(tree_pam6, names=spp, remove.cells = TRUE)


# Angiosperm species
angio_species <- tibble(species = tree_pam6$Spe) %>%
  mutate(species2 = gsub("_"," ", species)) %>%
  left_join(spname, by = c("species2" = "species")) %>%
  filter(class == "Magnoliopsida")

# only keep angiosperm species in the species presence/absence matrix
tree_pam6 <- lets.subsetPAM(tree_pam6, names = angio_species$species, remove.cells = TRUE)


## remove grid-cells with richness < 5
# the values of raster cells with richness <5 are set as NA 
tree_pam6$Richness_Raster[tree_pam6$Rich<5] <- NA 
table(!is.na(values(tree_pam6$Rich))) #TRUE:3098; FALSE:8734

## remove grid-cells that are not on mainland
land_bhm <- spTransform(land, behrmann)
tree_pam6$Richness_Raster <- mask(tree_pam6$Rich, land_bhm)  #based on mainland
table(!is.na(values(tree_pam6$Rich))) #FALSE:2319;TRUE:9513

# remove the cells from the presence/absence matrix and species with no occurrences
id_cell <- rownames(tree_pam6[[1]]) %in% which(!is.na(values(tree_pam6$Rich))) 
id_species <- colSums(tree_pam6[[1]][id_cell, -c(1:2)]) >0
tree_pam6$Presence_and_Absence_Matrix <- tree_pam6$Presence_and_Absence_Matrix[id_cell, c(1, 2, which(id_species) + 2)]

# update species names
tree_pam6$Species_name <- colnames(tree_pam6[[1]])[-c(1:2)]


# generate the polygons with grid-cells in land
raster_world <- tree_pam6$Rich
raster_world[] <- 1
grid_world <- rasterToPolygons(raster_world) 
colnames(grid_world@data)[1] <- "ID"
grid_world$ID <- 1:length(raster_world)
grid_land <- crop(grid_world, land_bhm)


# Species richness map
plot(tree_pam6, world=FALSE, axes=FALSE, box=FALSE)
plot(land_bhm, add=TRUE)


#Save data
save(tree_pam6, land_bhm, grid_land, file = "data/tree_pam/tree_pam6_final.RDATA")
