###################################################################################################
## Regress phylogenetic turnover (or nestedness) aganist functional turnover (or nestedness) and get the residuals.
## And then drow maps showing the residuals 

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "XUWUBING-PC" = "D:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

library(tidyverse)

source("Rcode/functions/mytmap.R")

load("models/beta_beta.deviation_envs.RDATA")
load("data/tree_pam/tree_pam6_final.RDATA")

# regress turnover and nestedness of functional beta-diversity agains respective components of phylogenetic beta-diversity
beta_env <- beta_env %>%
  mutate(resid.phylo.turn = resid(lm(phylo.beta.turn ~ func.beta.turn)),
         resid.phylo.nest = resid(lm(phylo.beta.nest ~ func.beta.nest))) %>%
  as.data.frame()

# the template of raster for mapping 
raster_cells <- tree_pam6$Richness_Raster
raster_cells[] <- NA

cor(beta_env[, c(29:30, 21:28)])


# Generate maps
# titles <- c(expression("A Residuals of FBD"[turnover]*" aganist PBD"[turnover]),
#            expression("B Residuals of FBD"[nestedness]*" aganist PBD"[nestedness]))
titles <- c("A Residuals of phylogenetic aganist functional turnover",
            "B Residuals of phylogenetic aganist functional nestedness")

resid.func.turn.map <- mytmap(myraster=raster_cells, myvalue=beta_env[,c(1, 29)], polygons=land_bhm, break.type="linear", mymidpoint=0,
                              digits=2,limit=0.01, style="cont", inner.margins=c(0.01,0,0.15,0), 
                              title=titles[1], title.size=0.7, title.position = c(0.10, 0.90),
                              legend.position=c(0.07, 0.01), legend.height=-0.6, legend.text.size=0.6)

resid.func.nest.map <- mytmap(myraster=raster_cells, myvalue=beta_env[,c(1, 30)], polygons=land_bhm, break.type="linear", mymidpoint=0,
                              digits=2,limit=0.01, style="cont", inner.margins=c(0.01,0,0.15,0), 
                              title=titles[2], title.size=0.7, title.position = c(0.10, 0.90),
                              legend.position=c(0.07, 0.01), legend.height=-0.6, legend.text.size=0.6)

resid.fbd.maps <- tmap_arrange(resid.func.turn.map, resid.func.nest.map, 
                               nrow=2, ncol=1, outer.margins=0)
tmap_save(resid.fbd.maps, file="results/Fig.S_residuals of PBD-turnover-nestedness.png",unit="mm",width=90, height=80)
tmap_save(resid.fbd.maps, file="results/Fig.S_residuals of PBD-turnover-nestedness.pdf",unit="mm",width=90, height=80)
