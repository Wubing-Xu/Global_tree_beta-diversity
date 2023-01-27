############################################################################################################################
### pre-statistical analyses: determine the predictor variables, whether and how transform response and predictor variables
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
                             mat = bio_1, map = bio_12, ts = bio_4, ps = bio_15, topo = topography, hmi = hmi)) %>%
  relocate(x, y, longitude, latitude, .after = gamma)

summary(beta_env)

# transform some variables
beta_env <- beta_env %>% 
  mutate(map.anomaly1 = abs(map.anomaly), # use absolute changes in precipitation indicating stability
         log_mat.anomaly = log10(mat.anomaly),
         log_map.anomaly = log10(map.anomaly1 + 1),
         log_ts = log10(ts),
         log_ps = log10(ps),
         log_map = log10(map),
         log_topo = log10(topo),
         sqrt_hmi = sqrt(hmi))

beta_env_scale <- beta_env %>% as.data.frame()
beta_env_scale[,9:36] <- scale(beta_env_scale[,9:36])


# Beta deviation from expected using local null model
betadev_env <- bind_cols(phylosor_devia[, 1:8], funcsor_devia[, 5:8]) %>% 
  left_join(beta_env %>% dplyr::select(ID, x:latitude, mat.anomaly:sqrt_hmi)) %>%
  relocate(x, y, longitude, latitude, .after = gamma)

betadev_env_scale <- betadev_env %>% as.data.frame()
betadev_env_scale[,9:32] <- scale(betadev_env_scale[,9:32])
  


#########
# evaluate whether response and predictor variables are needed to be transformed using simple and multiple linear regression

source("Rcode/functions/scatter.plot.R")

beta_env <- beta_env %>% as.data.frame()
betadev_env <- betadev_env %>% as.data.frame()

# A function to explore bivariate association to  check whether data-tranformation
Assump.Check <- function(x, y){
  par(mfrow = c(2,4))
  id <- !is.na(y)
  y <- y[id]; x<-x[id]
  hist(x)
  scatter.plot(x=x, y=y, xlab="x", ylab="y", text.pos="bottomright")
  lm1 <- lm(y~x)
  text(min(x), max(y), labels=c(paste("AIC = ", round(AIC(lm1),2), "\nR2=", round(summary(lm1)$r.squared,3))), pos=4)
  plot(lm1)
  y.resid <- resid(lm1)
  hist(y.resid)
}


#######
## histogram of response and predictor variables

# beta-diversity
par(mfrow=c(3, 4))
for (i in c(9:20)) hist(beta_env[,i], main=colnames(beta_env)[i], cex=0.8)

# deviation of phylogenetic and functional beta diversity
par(mfrow=c(2, 4))
for (i in c(9:16)) hist(betadev_env[,i], main=colnames(betadev_env)[i], cex=0.8)

# compare three transformation for beta-diversity:logit-, arcsin-root- and log- transformation
logit <- function(x) log(x/(1-x))
y1 <- beta_env[, 12]
y2 <- log(y1/(1-y1))
y3 <- asin(sqrt(y1))
y4 <-log(y1)
cor(y1, y2); cor(y1, y3); cor(y1, y4)
par(mfrow=c(2,2)); plot(y1,y2); plot(y1,y3); plot(y1,y4); plot(y2,y4)
# usually strongly correlated


# environmental variables
par(mfrow=c(3,3))
for (i in 21:29) hist(beta_env[,i], main=colnames(beta_env)[i], cex=0.8)

par(mfrow=c(3, 3))
for (i in c(30:36)) hist(beta_env[,i], main=colnames(beta_env)[i], cex=0.8)
# these variables are right-skewed distributed: mat_anomaly, map_anomaly1, map, topo, hmi


##########
## compare using simple linear regression
# check MAT
# turnover: log-log or level-level; nestedness: level-level; pnest: level-level
i = 23
j = 10
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i] + 7), y=beta_env[,j])
Assump.Check(x=log10(beta_env[,i] + 7), y=log10(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i] + 7), y=logit(beta_env[,j]))


# check MAP
# turnover: weak association; nestedness: level-log (logit, level); pnest: level-level (log, logit)
i = 24
j = 11
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i]), y=beta_env[,j])
Assump.Check(x=log10(beta_env[,i]), y=log10(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i]), y=logit(beta_env[,j]))


# check MAT anomaly
# turnover: level-log (or level); nestedness: level-level; pnest: level-level
i = 21
j = 12
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i]), y=beta_env[,j])
Assump.Check(x=log10(beta_env[,i]), y=log10(beta_env[,j]))
Assump.Check(x=log10(beta_env[,i]), y=logit(beta_env[,j]))


# check MAP anomaly
# turnover: absolute-log (or level); nestedness: absolute-level; pnest: absolute-level
i = 22 #22, 28
j = 10
Assump.Check(x = beta_env[,28], y = beta_env[,j])
Assump.Check(x=beta_env[,28], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,28], y=logit(beta_env[,j]))

Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])+1), y=beta_env[,j])
Assump.Check(x=log10(abs(beta_env[,i])+1), y=log10(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])+1), y=logit(beta_env[,j]))


# check temperature seasonality
# turnover: level-log (or level); nestedness: level-level; pnest: level-level
i = 25
j = 10
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=beta_env[,j])
Assump.Check(x=log10(abs(beta_env[,i])), y=log10(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=logit(beta_env[,j]))


# check precipitation seasonality
# turnover: log-level (or log); nestedness: weak association; pnest: weak association
i = 26
j = 11
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=beta_env[,j])
Assump.Check(x=log10(abs(beta_env[,i])), y=log10(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=logit(beta_env[,j]))


# check elevation range
# turnover: log-logit; nestedness: weak association; pnest: weak association
i = 27
j = 10
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=beta_env[,j])
Assump.Check(x=log10(abs(beta_env[,i])), y=log10(beta_env[,j]))
Assump.Check(x=log10(abs(beta_env[,i])), y=logit(beta_env[,j]))


# check Human Modification Index
# turnover: weak association; nestedness: weak association; pnest: weak association
i = 28
j = 12
Assump.Check(x = beta_env[,i], y = beta_env[,j])
Assump.Check(x=beta_env[,i], y=log10(beta_env[,j]))
Assump.Check(x=beta_env[,i], y=logit(beta_env[,j]))
Assump.Check(x=sqrt(beta_env[,i]), y=beta_env[,j])
Assump.Check(x=sqrt(beta_env[,i]), y=log10(beta_env[,j]))
Assump.Check(x=sqrt(beta_env[,i]), y=logit(beta_env[,j]))


##########
## correlations among variables
library(corrplot)

pre_cor <- cor(beta_env[, c(21:36)])

# MAT and temperature seasonality are strongly correlated (r = -0.91)
corrplot(pre_cor, method = 'number') 

plot(abs(beta_env[, c(21:22)]))
plot(beta_env$mat.anomaly, abs(beta_env$map.anomaly))
plot(beta_env$mat.anomaly, log(abs(beta_env$map.anomaly)+ 1))
plot(beta_env$mat.anomaly, abs(beta_env$map.anomaly)/beta_env$map)
plot(beta_env$mat.anomaly, beta_env$map.anomaly/beta_env$map)

cor(beta_env$mat.anomaly, log(abs(beta_env$map.anomaly)+ 1))


##########
## compare using multiple linear regression

## species turnover
lm.spe.turn <- lm(spe.beta.turn ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                    map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = beta_env)
summary(lm.spe.turn)

# compare subset models: avoid strongly correlated variables occurred in the same model
mm.lm.spe.turn <- dredge(lm.spe.turn, rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                           !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                           !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.spe.turn, delta < 10)
# note to choose: mat, log_map, mat.anm, map.anm1 (or map.anm), ts (missing), log_ps, top, sqrt_hmi


## species nestedness
lm.spe.nest <- lm(spe.beta.nest ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                    map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = beta_env)
summary(lm.spe.nest)

mm.lm.spe.nest <- dredge(lm.spe.nest, rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                           !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                           !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.spe.nest, delta < 10)
# note to choose: mat (missing), map, mat.anm, map.anm1, ts, ps (or log ps), top (or log_top), hmi (or sqrt_hmi)


## species beta
lm.spe.total <- lm(spe.beta.total ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                     map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = beta_env)
summary(lm.spe.total)

mm.lm.spe.total <- dredge(lm.spe.total, rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                            !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                            !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.spe.total, delta < 20)
# note to choose: mat, log_map, log_mat.anm, map.anm, ts (missing), log_ps, top, sqrt_hmi


## species proportion of nestedness
lm.spe.pnest <- lm(spe.pnest ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                     map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = beta_env)
summary(lm.spe.pnest)
mm.lm.spe.pnest <- dredge(lm.spe.pnest,rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                            !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                            !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.spe.pnest, delta < 10)
# note to choose: mat or ts, map, mat.anm, map.anm1, ps (or log_ps), log_top (or top), hmi (or sqrt_hmi)


## deviation of phylogenetic turnover 
lm.phylo.turn.dev <- lm(phylo.beta.turn  ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                    map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = betadev_env)
summary(lm.phylo.turn.dev)

# compare subset models: avoid strongly correlated variabels occurred in the same model
mm.lm.phylo.turn.dev <- dredge(lm.phylo.turn.dev, rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                                 !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                                 !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.phylo.turn.dev, delta < 10)
# note to choose: mat, map, mat.anm, log_map.anm (or map.anm), ts (missing), log_ps (or ps), log_top (or top), sqrt_hmi


## deviation of phylogenetic nestedness
lm.phylo.nest.dev <- lm(phylo.beta.nest  ~ mat.anomaly + log_mat.anomaly + map.anomaly + map.anomaly1 + log_map.anomaly + mat + 
                          map + log_map + ts + log_ts + ps + log_ps + topo + log_topo + hmi + sqrt_hmi, data = betadev_env)
summary(lm.phylo.nest.dev)

# compare subset models: avoid strongly correlated variabels occurred in the same model
mm.lm.phylo.nest.dev <- dredge(lm.phylo.nest.dev, rank="AIC", subset = !(mat.anomaly && log_mat.anomaly) & !(hmi && sqrt_hmi) & 
                                 !(map.anomaly && log_map.anomaly) & !(map.anomaly1 && log_map.anomaly) & !(map.anomaly && map.anomaly1) &
                                 !(map && log_map) & !(ts && log_ts) & !(ps && log_ps) & !(topo && log_topo) & !(mat && ts) & !(mat && log_ts))
subset(mm.lm.phylo.nest.dev, delta < 10)
# note to choose: mat, map, mat.anm, log_map.anm (map.anm1 or ,anm), log_ps, log_top, sqrt_hmi

## decision: although the best choice is not-consistent across response variables, the overall best choice is: mat, map, mat_anm, map_anm1, ts, ps, log_topo, sqrt_hmi




################################
## compare models with temperature seasonality and models without seasonality
 
## models with temperature seasonality 
nn <- dnearneigh(as.matrix(beta_env_scale[,7:8]), d1=0, d2=300, longlat=TRUE)
w <- nb2listw(nn, zero.policy=TRUE)
id1 <- c(21, 28, 23, 24, 25, 26, 34) # 21, 23, 28, 26, 29

# taxonomic total
mydata <- beta_env_scale[, c(9, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_spe.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_spe.total)

# taxonomic turnover
mydata <- beta_env_scale[, c(10, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_spe.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_spe.turn)

# taxonomic nestedness
mydata <- beta_env_scale[, c(11, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_spe.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_spe.nest)

# taxonomic proportion of nestedness
mydata <- beta_env_scale[, c(12, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_spe.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_spe.pnest)

# phylogenetic total
mydata <- beta_env_scale[, c(13, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.total)

# phylogenetic turnover
mydata <- beta_env_scale[, c(14, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.turn)

# phylogenetic nestedness
mydata <- beta_env_scale[, c(15, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.nest)

# phylogenetic proportion of nestedness
mydata <- beta_env_scale[, c(16, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.pnest)

# functional total
mydata <- beta_env_scale[, c(17, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.total)

# functional turnover
mydata <- beta_env_scale[, c(18, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.turn)

# functional nestedness
mydata <- beta_env_scale[, c(19, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.nest)

# functional proportion of nestedness
mydata <- beta_env_scale[, c(20, id1)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.pnest)


# for deviation of phylogenetic and functional beta
id2 <- c(17, 24, 19, 20, 21, 22, 30) # 21, 23, 28, 26, 29
# deviation of phylogenetic turnover
mydata <- betadev_env_scale[, c(10, id2)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.turn.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.turn.dev)

# deviation of phylogenetic nestedness
mydata <- betadev_env_scale[, c(11, id2)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_phylo.nest.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_phylo.nest.dev)

# deviation of functional turnover
mydata <- betadev_env_scale[, c(14, id2)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.turn.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.turn.dev)

# deviation of functional nestedness
mydata <- betadev_env_scale[, c(15, id2)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar_func.nest.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar_func.nest.dev)



###############
## models without temperature seasonality 

id3 <- c(21, 28, 23, 24, 26, 34) # 21, 23, 28, 26, 29
# taxonomic total
mydata <- beta_env_scale[, c(9, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_spe.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_spe.total)

# taxonomic turnover
mydata <- beta_env_scale[, c(10, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_spe.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_spe.turn)

# taxonomic nestedness
mydata <- beta_env_scale[, c(11, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_spe.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_spe.nest)

# taxonomic proportion of nestedness
mydata <- beta_env_scale[, c(12, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_spe.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_spe.pnest)

# phylogenetic total
mydata <- beta_env_scale[, c(13, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.total)

# phylogenetic turnover
mydata <- beta_env_scale[, c(14, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.turn)

# phylogenetic nestedness
mydata <- beta_env_scale[, c(15, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.nest)

# phylogenetic proportion of nestedness
mydata <- beta_env_scale[, c(16, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.pnest)

# functional total
mydata <- beta_env_scale[, c(17, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.total <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.total)

# functional turnover
mydata <- beta_env_scale[, c(18, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.turn <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.turn)

# functional nestedness
mydata <- beta_env_scale[, c(19, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.nest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.nest)

# functional proportion of nestedness
mydata <- beta_env_scale[, c(20, id3)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.pnest <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.pnest)


# for deviation of phylogenetic and functional beta
id4 <- c(17, 24, 19, 20, 22, 30) # 21, 23, 28, 26, 29
# deviation of phylogenetic turnover
mydata <- betadev_env_scale[, c(10, id4)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.turn.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.turn.dev)

# deviation of phylogenetic nestedness
mydata <- betadev_env_scale[, c(11, id4)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_phylo.nest.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_phylo.nest.dev)

# deviation of functional turnover
mydata <- betadev_env_scale[, c(14, id4)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.turn.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.turn.dev)

# deviation of functional nestedness
mydata <- betadev_env_scale[, c(15, id4)]
fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(colnames(mydata)[-1], collapse= "+")))
msar1_func.nest.dev <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
summary(msar1_func.nest.dev)



###############
## correlations between temperature anomaly and beta-diversity variables using modified t-test

for(i in 9:20) {
  mtest_res <- modified.ttest(x = beta_env_scale[, i], y = beta_env_scale[, "mat.anomaly"], coords = betadev_env_scale[, 7:8], nclass = 13) 
  print(mtest_res)
}

for(i in 9:16) {
  mtest_res <- modified.ttest(x = betadev_env_scale[, i], y = beta_env_scale[, "mat.anomaly"], coords = betadev_env_scale[, 7:8], nclass = 13) 
  print(mtest_res)
}
