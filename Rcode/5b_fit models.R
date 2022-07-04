############################################################################################################################
## Statistical analyses: both OLS and SAR models to analyze relationships between beta and beta-deviation and environmental variables
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
needed_libs <- c("tidyverse", "dplyr", "spdep", "spatialreg", "ncf", "MuMIn", "SpatialPack") 

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

load("intermediate_results/beta_observed_expected_deviation.RDATA")
load("intermediate_results/environments_final.RDATA")


# combine beta diversity and environmental variables
beta_env <- bind_cols(spesor[, 1:8], phylosor[, 5:8], funcsor[, 5:8]) %>%
  # add environmental variables
  inner_join(env200_mean_nn24 %>% 
               dplyr::select(ID, x, y, longitude, latitude, mat.anomaly, map.anomaly, 
                             mat = bio_1, map = bio_12, ts = bio_4, ps = bio_15, topo = topography)) %>%
  relocate(x, y, longitude, latitude, .after = gamma)

summary(beta_env)

# transform some variables
beta_env <- beta_env %>%
  mutate(map.anomaly = abs(map.anomaly), # use absolute precipitation changes indicating stability
         log_topo = log10(topo))

beta_env_scale <- beta_env %>% as.data.frame()
beta_env_scale[,9:28] <- scale(beta_env_scale[,9:28])


# Beta deviation from expected using local null model
betadev_env <- bind_cols(phylosor_devia[, 1:8], funcsor_devia[, 5:8]) %>% 
  left_join(beta_env %>% dplyr::select(ID, x:latitude, mat.anomaly:log_topo)) %>%
  relocate(x, y, longitude, latitude, .after = gamma)

betadev_env_scale <- betadev_env %>% as.data.frame()
betadev_env_scale[,9:24] <- scale(betadev_env_scale[,9:24])
  
save(beta_env, betadev_env, file="models/beta_beta.deviation_envs.RDATA")


###################
## Use multiple OLS to calculate relative importance of predictors
# Observed beta
mm_ols_spe.total <- OLScoef(mydata = beta_env_scale[, c(9, 21:26, 28)]) 
mm_ols_spe.turn <- OLScoef(mydata = beta_env_scale[, c(10, 21:26, 28)]) 
mm_ols_spe.nest <- OLScoef(mydata = beta_env_scale[, c(11, 21:26, 28)]) 
mm_ols_spe.pnest <- OLScoef(mydata = beta_env_scale[, c(12, 21:26, 28)]) 

mm_ols_phylo.total <- OLScoef(mydata = beta_env_scale[, c(13, 21:26, 28)]) 
mm_ols_phylo.turn <- OLScoef(mydata = beta_env_scale[, c(14, 21:26, 28)]) 
mm_ols_phylo.nest <- OLScoef(mydata = beta_env_scale[, c(15, 21:26, 28)]) 
mm_ols_phylo.pnest <- OLScoef(mydata = beta_env_scale[, c(16, 21:26, 28)]) 

mm_ols_func.total <- OLScoef(mydata = beta_env_scale[, c(17, 21:26, 28)]) 
mm_ols_func.turn <- OLScoef(mydata = beta_env_scale[, c(18, 21:26, 28)]) 
mm_ols_func.nest <- OLScoef(mydata = beta_env_scale[, c(19, 21:26, 28)]) 
mm_ols_func.pnest <- OLScoef(mydata = beta_env_scale[, c(20, 21:26, 28)]) 

# Beta deviation from expected using local null model
mm_ols_phylo.dev.total <- OLScoef(mydata = betadev_env_scale[, c(9, 17:22, 24)]) 
mm_ols_phylo.dev.turn <- OLScoef(mydata = betadev_env_scale[, c(10, 17:22, 24)]) 
mm_ols_phylo.dev.nest <- OLScoef(mydata = betadev_env_scale[, c(11, 17:22, 24)]) 
mm_ols_phylo.dev.pnest <- OLScoef(mydata = betadev_env_scale[, c(12, 17:22, 24)]) 
mm_ols_func.dev.total <- OLScoef(mydata = betadev_env_scale[, c(13, 17:22, 24)]) 
mm_ols_func.dev.turn <- OLScoef(mydata = betadev_env_scale[, c(14, 17:22, 24)]) 
mm_ols_func.dev.nest <- OLScoef(mydata = betadev_env_scale[, c(15, 17:22, 24)]) 
mm_ols_func.dev.pnest <- OLScoef(mydata = betadev_env_scale[, c(16, 17:22, 24)]) 

save(mm_ols_spe.total, mm_ols_spe.turn, mm_ols_spe.nest, mm_ols_spe.pnest,
     mm_ols_phylo.total, mm_ols_phylo.turn, mm_ols_phylo.nest, mm_ols_phylo.pnest,
     mm_ols_func.total, mm_ols_func.turn, mm_ols_func.nest, mm_ols_func.pnest,
     mm_ols_phylo.dev.total, mm_ols_phylo.dev.turn, mm_ols_phylo.dev.nest, mm_ols_phylo.dev.pnest,
     mm_ols_func.dev.total, mm_ols_func.dev.turn, mm_ols_func.dev.nest, mm_ols_func.dev.pnest,
     file="models/MM_OLS_beta_beta.deviation_envs.RDATA")



#########
##SAR analysis

## Determine the appropriate distance to do SAR analysis
# co.increment is for the correlogram. 250 km is the maximum value of distance between a focal cell and its nearest neighbors 
# observed beta
distSAR_spe.total <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,9], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_spe.turn <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,10], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_spe.nest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,11], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_spe.pnest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,12], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km

distSAR_phylo.total <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,13], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_phylo.turn <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,14], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_phylo.nest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,15], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_phylo.pnest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,16], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km

distSAR_func.total <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,17], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.turn <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,18], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.nest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,19], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.pnest <- distSAR(coords=beta_env_scale[,7:8], y=beta_env_scale[,20], env=beta_env_scale[, c(21:26, 28)], start=200, end=1000, co.increment=250) #300 km

# beta deviation
distSAR_phylo.dev.total <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,9], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km
distSAR_phylo.dev.turn <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,10], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km
distSAR_phylo.dev.nest <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,11], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 or 200 km
distSAR_phylo.dev.pnest <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,12], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km

distSAR_func.dev.total <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,13], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.dev.turn <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,14], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.dev.nest <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,15], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km
distSAR_func.dev.pnest <- distSAR(coords=betadev_env_scale[,7:8], y=betadev_env_scale[,16], env=betadev_env_scale[,c(17:22, 24)], start=200, end=1000, co.increment=250) #300 km

save(distSAR_spe.total, distSAR_spe.turn, distSAR_spe.nest, distSAR_spe.pnest,
     distSAR_phylo.total, distSAR_phylo.turn, distSAR_phylo.nest, distSAR_phylo.pnest,
     distSAR_func.total, distSAR_func.turn, distSAR_func.nest, distSAR_func.pnest,
     distSAR_phylo.dev.total, distSAR_phylo.dev.turn, distSAR_phylo.dev.nest, distSAR_phylo.dev.pnest,
     distSAR_func.dev.total, distSAR_func.dev.turn, distSAR_func.dev.nest, distSAR_func.dev.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs_best.dist.RDATA")


## Use multiple SAR to calculate relative importance of predictors
# observed beta
mm_sar_spe.total <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(9, 21:26, 28)], sardist=300) 
mm_sar_spe.turn <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(10, 21:26, 28)], sardist=300) 
mm_sar_spe.nest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(11, 21:26, 28)], sardist=300) 
mm_sar_spe.pnest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(12, 21:26, 28)], sardist=300) 

save(mm_sar_spe.total, mm_sar_spe.turn, mm_sar_spe.nest, mm_sar_spe.pnest,
    file="models/MM_SAR_beta_beta.deviation_envs_1.RDATA")

mm_sar_phylo.total <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(13, 21:26, 28)], sardist=300) 
mm_sar_phylo.turn <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(14, 21:26, 28)], sardist=300) 
mm_sar_phylo.nest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(15, 21:26, 28)], sardist=300) 
mm_sar_phylo.pnest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(16, 21:26, 28)], sardist=300) 

save(mm_sar_phylo.total, mm_sar_phylo.turn, mm_sar_phylo.nest, mm_sar_phylo.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs_2.RDATA")

mm_sar_func.total <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(17, 21:26, 28)], sardist=300) 
mm_sar_func.turn <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(18, 21:26, 28)], sardist=300) 
mm_sar_func.nest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(19, 21:26, 28)], sardist=300) 
mm_sar_func.pnest <- MSARcoef(coords=beta_env_scale[,7:8], mydata=beta_env_scale[, c(20, 21:26, 28)], sardist=300) 

save(mm_sar_func.total, mm_sar_func.turn, mm_sar_func.nest, mm_sar_func.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs_3.RDATA")

# beta deviation from expected using local null model
mm_sar_phylo.dev.total <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(9, 17:22, 24)], sardist=300) 
mm_sar_phylo.dev.turn <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(10, 17:22, 24)], sardist=300) 
mm_sar_phylo.dev.nest <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(11, 17:22, 24)], sardist=300) 
mm_sar_phylo.dev.pnest <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(12, 17:22, 24)], sardist=300) 

save(mm_sar_phylo.dev.total, mm_sar_phylo.dev.turn, mm_sar_phylo.dev.nest, mm_sar_phylo.dev.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs_4.RDATA")

mm_sar_func.dev.total <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(13, 17:22, 24)], sardist=300)
mm_sar_func.dev.turn <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(14, 17:22, 24)], sardist=300)
mm_sar_func.dev.nest <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(15, 17:22, 24)], sardist=300)
mm_sar_func.dev.pnest <- MSARcoef(coords=betadev_env_scale[,7:8], mydata=betadev_env_scale[, c(16, 17:22, 24)], sardist=300)

save(mm_sar_func.dev.total, mm_sar_func.dev.turn, mm_sar_func.dev.nest, mm_sar_func.dev.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs_5.RDATA")


save(mm_sar_spe.total, mm_sar_spe.turn, mm_sar_spe.nest, mm_sar_spe.pnest,
     mm_sar_phylo.total, mm_sar_phylo.turn, mm_sar_phylo.nest, mm_sar_phylo.pnest,
     mm_sar_func.total, mm_sar_func.turn, mm_sar_func.nest, mm_sar_func.pnest,
     mm_sar_phylo.dev.total, mm_sar_phylo.dev.turn, mm_sar_phylo.dev.nest, mm_sar_phylo.dev.pnest,
     mm_sar_func.dev.total, mm_sar_func.dev.turn, mm_sar_func.dev.nest, mm_sar_func.dev.pnest,
     file="models/MM_SAR_beta_beta.deviation_envs.RDATA")
