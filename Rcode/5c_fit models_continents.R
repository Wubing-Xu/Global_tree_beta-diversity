############################################################################################################################
## use OLS and SAR models to analyze relationships between beta and environmental variables within each of six continents
############################################################################################################################

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "XUWUBING-PC" = "D:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse", "dplyr", "spdep", "spatialreg", "ncf", "MuMIn", "SpatialPack", "rworldmap") 

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)
options(na.action="na.fail")

# input self-defined functions to run OLS and SAR
source("Rcode/functions/distSAR.R")
source("Rcode/functions/Multiple_OLS_SAR.R")

# load the beta and environmental data
load("models/beta_beta.deviation_envs.RDATA")


# get the global map with the information of continents 
country <- getMap(resolution='low')
plot(country, col = as.factor(country@data$REGION))

# set the Russia is in Asia
country@data$REGION <- as.character(country@data$REGION)
country@data$REGION[country@data$NAME == "Russia" & !is.na(country@data$NAME)] = "Asia"
table(country$REGION)
plot(country, col = as.factor(country@data$REGION))


# get the continents that each grid cell located
grid_points <-  SpatialPoints(data.frame(beta_env$longitude, beta_env$latitude), proj4string=CRS(proj4string(country)))  

beta_env <- beta_env %>% 
  mutate(continent = over(grid_points, country)$REGION)

table(beta_env$continent, useNA = "always")

# 6 grid cells match no continents; add them manually
beta_env %>% filter(is.na(continent)) %>% dplyr::select(ID, longitude, latitude, continent)
# ID_manual <- data.frame("ID" = c(1137, 1266, 2373, 2548, 2549, 4141, 5978),
#                        continent = c("Europe", "North America", "Europe", "Asia", "Asia", "Asia", "South America"))
ID_manual <- tibble("ID" = c(1093, 1137, 1266, 2373, 2548, 4141),
                    continent = c("North America","Europe", "North America", "Europe", "Asia", "Asia"))

id <- match(ID_manual$ID, beta_env$ID)
beta_env$continent[id] <- ID_manual$continent

# add the continents for "betadev_env"
betadev_env <- betadev_env %>% 
  left_join(beta_env %>% dplyr::select(ID, continent))

table(betadev_env$continent, useNA = "always")


# z-transform response and predictor variables, which will be used for models
beta_env_scale <- beta_env %>%
  group_by(continent) %>%
  mutate(across(spe.beta.total:sqrt_hmi, scale)) %>%
  ungroup() %>%
  as.data.frame()

betadev_env_scale <- betadev_env %>%
  group_by(continent) %>%
  mutate(across(phylo.beta.total:sqrt_hmi, scale)) %>%
  ungroup() %>%
  as.data.frame()


#######################
## fit models
## Use multiple OLS to calculate relative importance of predictors
mm_ols_region_spe.total <- list()
mm_ols_region_spe.turn <- list()
mm_ols_region_spe.nest <- list()
mm_ols_region_spe.pnest <- list()
mm_ols_region_phylo.total <- list()
mm_ols_region_phylo.turn <- list()
mm_ols_region_phylo.nest <- list()
mm_ols_region_phylo.pnest <- list()
mm_ols_region_func.total <- list()
mm_ols_region_func.turn <- list()
mm_ols_region_func.nest <- list()
mm_ols_region_func.pnest <- list()
mm_ols_region_phylo.dev.turn <- list()
mm_ols_region_phylo.dev.nest <- list()
mm_ols_region_func.dev.turn <- list()
mm_ols_region_func.dev.nest  <- list()

continents <- unique(beta_env$continent)

for(i in 1:length(continents)){
  continent <- continents[i]
  beta_env_sub <- beta_env_scale[beta_env_scale$continent == continent, ]
  betadev_env_sub <- betadev_env_scale[betadev_env_scale$continent == continent, ]
  
  # Observed beta
  mm_ols_region_spe.total[[i]] <- OLScoef(mydata = beta_env_sub[, c(9, 21:26, 29, 30)]) 
  mm_ols_region_spe.turn[[i]] <- OLScoef(mydata = beta_env_sub[, c(10, 21:26, 29, 30)]) 
  mm_ols_region_spe.nest[[i]] <- OLScoef(mydata = beta_env_sub[, c(11, 21:26, 29, 30)]) 
  mm_ols_region_spe.pnest[[i]] <- OLScoef(mydata = beta_env_sub[, c(12, 21:26, 29, 30)]) 
  
  mm_ols_region_phylo.total[[i]] <- OLScoef(mydata = beta_env_sub[, c(13, 21:26, 29, 30)]) 
  mm_ols_region_phylo.turn[[i]] <- OLScoef(mydata = beta_env_sub[, c(14, 21:26, 29, 30)]) 
  mm_ols_region_phylo.nest[[i]] <- OLScoef(mydata = beta_env_sub[, c(15, 21:26, 29, 30)]) 
  mm_ols_region_phylo.pnest[[i]] <- OLScoef(mydata = beta_env_sub[, c(16, 21:26, 29, 30)]) 
  
  mm_ols_region_func.total[[i]] <- OLScoef(mydata = beta_env_sub[, c(17, 21:26, 29, 30)]) 
  mm_ols_region_func.turn[[i]] <- OLScoef(mydata = beta_env_sub[, c(18, 21:26, 29, 30)]) 
  mm_ols_region_func.nest[[i]] <- OLScoef(mydata = beta_env_sub[, c(19, 21:26, 29, 30)]) 
  mm_ols_region_func.pnest[[i]] <- OLScoef(mydata = beta_env_sub[, c(20, 21:26, 29, 30)]) 
  
  # Beta deviation from expected using local null model
  mm_ols_region_phylo.dev.turn[[i]] <- OLScoef(mydata = betadev_env_sub[, c(10, 17:22, 25, 26)]) 
  mm_ols_region_phylo.dev.nest[[i]] <- OLScoef(mydata = betadev_env_sub[, c(11, 17:22, 25, 26)]) 
  mm_ols_region_func.dev.turn[[i]] <- OLScoef(mydata = betadev_env_sub[, c(14, 17:22, 25, 26)]) 
  mm_ols_region_func.dev.nest[[i]] <- OLScoef(mydata = betadev_env_sub[, c(15, 17:22, 25, 26)]) 
  
  print(i)
}

names(mm_ols_region_spe.total) <- continents; names(mm_ols_region_spe.turn) <- continents
names(mm_ols_region_spe.nest) <- continents; names(mm_ols_region_spe.pnest) <- continents

names(mm_ols_region_phylo.total) <- continents; names(mm_ols_region_phylo.turn) <- continents
names(mm_ols_region_phylo.nest) <- continents; names(mm_ols_region_phylo.pnest) <- continents

names(mm_ols_region_func.total) <- continents; names(mm_ols_region_func.turn) <- continents
names(mm_ols_region_func.nest) <- continents; names(mm_ols_region_func.pnest) <- continents

names(mm_ols_region_phylo.dev.turn) <- continents; names(mm_ols_region_phylo.dev.nest) <- continents
names(mm_ols_region_func.dev.turn) <- continents; names(mm_ols_region_func.dev.nest) <- continents

save(mm_ols_region_spe.total, mm_ols_region_spe.turn, mm_ols_region_spe.nest, mm_ols_region_spe.pnest,
     mm_ols_region_phylo.total, mm_ols_region_phylo.turn, mm_ols_region_phylo.nest, mm_ols_region_phylo.pnest,
     mm_ols_region_func.total, mm_ols_region_func.turn, mm_ols_region_func.nest, mm_ols_region_func.pnest,
     mm_ols_region_phylo.dev.turn, mm_ols_region_phylo.dev.nest, mm_ols_region_func.dev.turn, mm_ols_region_func.dev.nest,
     file="models/MM_OLS_region_beta_beta.deviation_envs.RDATA")


#######################
## fit models
## Use multiple SAR to calculate relative importance of predictors

mm_sar_region_spe.total <- list()
mm_sar_region_spe.turn <- list()
mm_sar_region_spe.nest <- list()
mm_sar_region_spe.pnest <- list()
mm_sar_region_phylo.total <- list()
mm_sar_region_phylo.turn <- list()
mm_sar_region_phylo.nest <- list()
mm_sar_region_phylo.pnest <- list()
mm_sar_region_func.total <- list()
mm_sar_region_func.turn <- list()
mm_sar_region_func.nest <- list()
mm_sar_region_func.pnest <- list()
mm_sar_region_phylo.dev.turn <- list()
mm_sar_region_phylo.dev.nest <- list()
mm_sar_region_func.dev.turn <- list()
mm_sar_region_func.dev.nest  <- list()

continents <- unique(beta_env$continent)

for(i in 1:length(continents)){
  continent <- continents[i]
  beta_env_sub <- beta_env_scale[beta_env_scale$continent == continent, ]
  betadev_env_sub <- betadev_env_scale[betadev_env_scale$continent == continent, ]
  
  mm_sar_region_spe.total[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(9, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_spe.turn[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(10, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_spe.nest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(11, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_spe.pnest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(12, 21:26, 29, 30)], sardist=300) 
  
  print(i)
}

names(mm_sar_region_spe.total) <- continents; names(mm_sar_region_spe.turn) <- continents
names(mm_sar_region_spe.nest) <- continents; names(mm_sar_region_spe.pnest) <- continents

save(mm_sar_region_spe.total, mm_sar_region_spe.turn, mm_sar_region_spe.nest, mm_sar_region_spe.pnest, 
     file="models/MM_SAR_beta_region_beta.deviation_envs_1.RDATA" )


for(i in 1:length(continents)){
  continent <- continents[i]
  beta_env_sub <- beta_env_scale[beta_env_scale$continent == continent, ]
  betadev_env_sub <- betadev_env_scale[betadev_env_scale$continent == continent, ]
  
  mm_sar_region_phylo.total[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(13, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_phylo.turn[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(14, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_phylo.nest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(15, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_phylo.pnest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(16, 21:26, 29, 30)], sardist=300) 
  
  print(i)
}

names(mm_sar_region_phylo.total) <- continents; names(mm_sar_region_phylo.turn) <- continents
names(mm_sar_region_phylo.nest) <- continents; names(mm_sar_region_phylo.pnest) <- continents

save(mm_sar_region_phylo.total, mm_sar_region_phylo.turn, mm_sar_region_phylo.nest, mm_sar_region_phylo.pnest, 
     file="models/MM_SAR_beta_region_beta.deviation_envs_2.RDATA" )


for(i in 1:length(continents)){
  continent <- continents[i]
  beta_env_sub <- beta_env_scale[beta_env_scale$continent == continent, ]
  betadev_env_sub <- betadev_env_scale[betadev_env_scale$continent == continent, ]
  
  mm_sar_region_func.total[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(17, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_func.turn[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(18, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_func.nest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(19, 21:26, 29, 30)], sardist=300) 
  mm_sar_region_func.pnest[[i]] <- MSARcoef(coords=beta_env_sub[,7:8], mydata=beta_env_sub[, c(20, 21:26, 29, 30)], sardist=300) 
  
  print(i)
}

names(mm_sar_region_func.total) <- continents; names(mm_sar_region_func.turn) <- continents
names(mm_sar_region_func.nest) <- continents; names(mm_sar_region_func.pnest) <- continents

save(mm_sar_region_func.total, mm_sar_region_func.turn, mm_sar_region_func.nest, mm_sar_region_func.pnest, 
     file="models/MM_SAR_beta_region_beta.deviation_envs_3.RDATA" )


for(i in 1:length(continents)){
  continent <- continents[i]
  beta_env_sub <- beta_env_scale[beta_env_scale$continent == continent, ]
  betadev_env_sub <- betadev_env_scale[betadev_env_scale$continent == continent, ]
  
  # Beta deviation from expected using local null model
  mm_sar_region_phylo.dev.turn[[i]] <- MSARcoef(coords=betadev_env_sub[,7:8], mydata=betadev_env_sub[, c(10, 17:22, 25, 26)], sardist=300) 
  mm_sar_region_phylo.dev.nest[[i]] <- MSARcoef(coords=betadev_env_sub[,7:8], mydata=betadev_env_sub[, c(11, 17:22, 25, 26)], sardist=300) 
  mm_sar_region_func.dev.turn[[i]] <- MSARcoef(coords=betadev_env_sub[,7:8], mydata=betadev_env_sub[, c(14, 17:22, 25, 26)], sardist=300) 
  mm_sar_region_func.dev.nest[[i]] <- MSARcoef(coords=betadev_env_sub[,7:8], mydata=betadev_env_sub[, c(15, 17:22, 25, 26)], sardist=300) 
  print(i)
}

names(mm_sar_region_phylo.dev.turn) <- continents; names(mm_sar_region_phylo.dev.nest) <- continents
names(mm_sar_region_func.dev.turn) <- continents; names(mm_sar_region_func.dev.nest) <- continents

save(mm_sar_region_phylo.dev.turn, mm_sar_region_phylo.dev.nest, mm_sar_region_func.dev.turn, mm_sar_region_func.dev.nest, 
     file="models/MM_SAR_beta_region_beta.deviation_envs_4.RDATA" )


save(mm_sar_region_spe.total, mm_sar_region_spe.turn, mm_sar_region_spe.nest, mm_sar_region_spe.pnest,
     mm_sar_region_phylo.total, mm_sar_region_phylo.turn, mm_sar_region_phylo.nest, mm_sar_region_phylo.pnest,
     mm_sar_region_func.total, mm_sar_region_func.turn, mm_sar_region_func.nest, mm_sar_region_func.pnest,
     mm_sar_region_phylo.dev.turn, mm_sar_region_phylo.dev.nest, mm_sar_region_func.dev.turn, mm_sar_region_func.dev.nest,
     beta_env, betadev_env, country, 
     file="models/MM_SAR_region_beta_beta.deviation_envs.RDATA")
