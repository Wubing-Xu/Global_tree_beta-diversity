################################################################################################
## Extract the environmental conditions for all grid cells in resolution of 200km;
## and then calculate the the mean of environmental conditions across 3*3 or 5*5 grid cells
################################################################################################

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","letsR", "raster", "spdep", "sp")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


## get spatial grid cells
load("data/tree_pam/tree_pam6_final.RDATA")

## load environmental variables
# elevation
elev_dir <-paste("data/environment_rasters/elevation/wc2.1_30s_elev.tif")
elev <- stack(elev_dir)

# current climates
bioc_dir <-paste("data/environment_rasters/current_climate/wc2.1_5m_bio/wc2.1_5m_bio_", 1:19, ".tif", sep="")
bioc <- stack(bioc_dir)

# LGM climates
lgmc_dir <-paste("data/environment_rasters/LGM_climate/chelsa_LGM_v1_2B_r5m/5min/bio_", c(1,12), ".tif", sep="")
lgmc <- stack(lgmc_dir)


# A function to extract the value of environments for each mypolygon 
extract_env <- function(env, mypolygon, res, fun=mean, weights=FALSE){
  CRS_mypolygon <- projection(mypolygon)

  env <- projectRaster(env, crs=CRS_mypolygon, res=res)
  env_mypolygon_value <- extract(env, mypolygon, fun=fun, weights=weights, na.rm=TRUE, df=TRUE) 
  return(env_mypolygon_value)
}

# Calculate current and LGM climates for each grid cell
bioc_grid <- extract_env(env=bioc, mypolygon=grid_land, res=10)
bioc_grid[, c(5)] <- bioc_grid[, 5]/100 #the raw unit is standard deviation*100
bioc_grid[, 1] <- grid_land@data[, 1]

lgmc_grid <- extract_env(env=lgmc, mypolygon=grid_land, res=10)
lgmc_grid[, c(2)] <- lgmc_grid[, 2]/10 #the raw unit is 1/10 degree
lgmc_grid[, 1] <- grid_land@data[, 1]

## Calculate mean and range of elevation
get_range <- function(x, na.rm=TRUE) {
  range = max(x, na.rm=na.rm) - min(x,na.rm=TRUE)
  return(range)
}
# elevational range
elev_range_grid200 <- extract_env(env=elev, mypolygon=grid_land, res=1, fun=get_range)
colnames(elev_range_grid200)[2] <- "topography"
elev_range_grid200[, 1] <- grid_land@data[, 1]

# mean elevatioins
elev_5m <- aggregate(elev, 10)
elev_mean_grid200 <- extract_env(env=elev_5m, mypolygon=grid_land, res=10, fun=mean)
colnames(elev_mean_grid200)[2] <- "elevation"
elev_mean_grid200[, 1] <- grid_land@data[, 1]


save(bioc_grid, lgmc_grid, elev_range_grid200, elev_mean_grid200,
     file="intermediate_results/environments_allCells.RDATA")
load("intermediate_results/environments_allCells.RDATA")



#################
## Assemble environmental variables

# the projected coordinates 
xy <- coordinates(tree_pam6[[2]])[grid_land@data[, 1], ]
colnames(xy) <- c("x", "y")

# get longitude and latitude
cell_points <- as.data.frame(xy)
coordinates(cell_points) <- c("x","y")
projection(cell_points) <- projection(tree_pam6[[2]])
cell_points_longlat <- spTransform(cell_points, CRS("+proj=longlat +datum=WGS84"))
long_lat <- coordinates(cell_points_longlat)
colnames(long_lat) <- c("longitude", "latitude")

# temperature and precipitation anomaly since the LGM
lgmcc_grid <- data.frame(ID = grid_land@data[, 1], 
                         mat.anomaly = bioc_grid[, 2] - lgmc_grid[, 2],
                         map.anomaly = bioc_grid[, 13] - lgmc_grid[, 3])

# change names of bioclimatic variables
colnames(bioc_grid)[-1] <- paste0("bio_", 1:19)

# combine environmental variables
env200 <- data.frame(ID = grid_land@data[, 1], 
                     xy, long_lat, bioc_grid[, -1], lgmcc_grid[, -1],
                     elevation = elev_mean_grid200[,-1],
                     topography = elev_range_grid200[,-1]) %>%
  as_tibble() %>%
  # remove data of grid-cells with small part in the land
  mutate(land_area = rgeos::gArea(grid_land, byid = TRUE)) %>%
  filter(land_area >= 4000) %>%
  dplyr::select(-land_area)



##############################################
##Calculate the mean of environmmental conditions of focal cells and their eight neighboring cells

# the environment subset of cells with tree distributions that were used to calculate beta
tree_cells <- which(!is.na(tree_pam6$Richness_Raster[]))
env200_treecell <- env200 %>% 
  filter(ID %in% tree_cells)

# define the 8 nearest neighboring cells for each focal cells, and include including itself
nb8 <- dnearneigh(x = as.matrix(env200_treecell[,2:3]), d1 = 0, d2 = 300, longlat = FALSE)
nb8mat <- nb2mat(neighbours = nb8, style = "B", zero.policy = TRUE)
diag(nb8mat) <- 1

# define the 24 nearest neighboring cells for each focal cells, and includeincluding itself
nb24 <- dnearneigh(x = as.matrix(env200_treecell[,2:3]), d1 = 0, d2 = 570, longlat = FALSE)
nb24mat <- nb2mat(neighbours = nb24, style = "B", zero.policy = TRUE)
diag(nb24mat) <- 1

# define the 24 nearest neighboring cells for each focal cells use all grid cells (not just cells that have tree observations)
nb24_all <- dnearneigh(x = as.matrix(env200[,2:3]), d1 = 0, d2 = 570, longlat = FALSE)
nb24mat_all <- nb2mat(neighbours = nb24_all, style = "B", zero.policy = TRUE)
diag(nb24mat_all) <- 1


# A function to calculate environmental conditions of focal and neighboring cells based on a function
get_env_nb_summ <- function(envdata = envdata, nbmat = NA, fun){
  envdata.res <- envdata # a new data frame to store results
  envdata.res[, -1] <- NA
  
  get_env_nb_cell <- function(x){
    id <- which(x == 1)
    if(ncol(envdata) == 2) env.nb.cell <- apply(as.data.frame(envdata[id,-1]), 2 ,fun, na.rm=TRUE)
    if(ncol(envdata) > 2)  env.nb.cell <- apply(envdata[id,-1], 2, fun, na.rm=TRUE)
    return(env.nb.cell)
  }
  env.nb.cells <- apply(nbmat, 1, get_env_nb_cell)
  if(ncol(envdata)>2) env.nb.cells <- t(env.nb.cells)
  
  envdata.res[, -1] <- env.nb.cells
  return(envdata.res)
}

# mean environmental conditions
env200_mean_nn8 <- env200_treecell
env200_mean_nn8[,c(1, 6:28)] <- get_env_nb_summ(envdata = env200_treecell[,c(1, 6:28)], nbmat = nb8mat, fun=mean)

env200_mean_nn24 <- env200_treecell
env200_mean_nn24[,c(1, 6:28)] <- get_env_nb_summ(envdata = env200_treecell[,c(1, 6:28)], nbmat = nb24mat, fun=mean)

env200_mean_nn24_all <- env200
env200_mean_nn24_all[,c(1, 6:28)] <- get_env_nb_summ(envdata = env200[,c(1, 6:28)], nbmat = nb24mat_all, fun=mean)

# save all calculated output
save(env200, env200_mean_nn8, env200_mean_nn24, env200_mean_nn24_all, file = "intermediate_results/environments_final.RDATA")
