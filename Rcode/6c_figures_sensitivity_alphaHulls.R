###################################################################################################
# draw figures of sensitivity analyses that range maps based on different alpha values (2, 4, 10) were used in calculating beta-diversity
# draw maps for results using alpha of 2, 4, and; and a correlation plot showing correlations of beta based on alpha 6 vervus 2, 4, and 10

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "XUWUBING-PC" = "D:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)

library(tidyverse)

# load packages
needed_libs <- c("tidyverse","ggplot2", "ggridges", "cowplot", "ggstance")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


source("Rcode/functions/mytmap.R")

load("intermediate_results/observed_beta_sorensen.RDATA")
load("data/tree_pam/tree_pam6_final.RDATA")
load("data/tree_pam/tree_pam2_final.RDATA")


## organize beta diversity variable using ranges that estimated with alpha hulls using different values of alphas
# alpha hull 2
beta_alpha2 <- bind_cols(
  spe.sor_pam2[, c(1, 7, 5, 6, 12)] %>% 
  as_tibble() %>%
  rename(spe_total_alpha2= beta.SOR, spe_turn_alpha2 = beta.SIM, spe_nest_alpha2 = beta.SNE, spe_pnest_alpha2 = pnest),
  phylo.sor_pam2[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_alpha2= phylo.beta.SOR, phylo_turn_alpha2 = phylo.beta.SIM, phylo_nest_alpha2 = phylo.beta.SNE, phylo_pnest_alpha2 = pnest),
  func.sor_pam2[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_alpha2= phylo.beta.SOR, func_turn_alpha2 = phylo.beta.SIM, func_nest_alpha2 = phylo.beta.SNE, func_pnest_alpha2 = pnest))

# alpha hull 4
beta_alpha4 <- bind_cols(
  spe.sor_pam4[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(spe_total_alpha4= beta.SOR, spe_turn_alpha4 = beta.SIM, spe_nest_alpha4 = beta.SNE, spe_pnest_alpha4 = pnest),
  phylo.sor_pam4[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_alpha4= phylo.beta.SOR, phylo_turn_alpha4 = phylo.beta.SIM, phylo_nest_alpha4 = phylo.beta.SNE, phylo_pnest_alpha4 = pnest),
  func.sor_pam4[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_alpha4= phylo.beta.SOR, func_turn_alpha4 = phylo.beta.SIM, func_nest_alpha4 = phylo.beta.SNE, func_pnest_alpha4 = pnest))

# alpha hull 6
beta_alpha6 <- bind_cols(
  spe.sor[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(spe_total_alpha6= beta.SOR, spe_turn_alpha6 = beta.SIM, spe_nest_alpha6 = beta.SNE, spe_pnest_alpha6 = pnest),
  phylo.sor[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_alpha6= phylo.beta.SOR, phylo_turn_alpha6 = phylo.beta.SIM, phylo_nest_alpha6 = phylo.beta.SNE, phylo_pnest_alpha6 = pnest),
  func.sor[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_alpha6= phylo.beta.SOR, func_turn_alpha6 = phylo.beta.SIM, func_nest_alpha6 = phylo.beta.SNE, func_pnest_alpha6 = pnest))

# alpha hull 10
beta_alpha10 <- bind_cols(
  spe.sor_pam10[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(spe_total_alpha10= beta.SOR, spe_turn_alpha10 = beta.SIM, spe_nest_alpha10 = beta.SNE, spe_pnest_alpha10 = pnest),
  phylo.sor_pam10[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_alpha10= phylo.beta.SOR, phylo_turn_alpha10 = phylo.beta.SIM, phylo_nest_alpha10 = phylo.beta.SNE, phylo_pnest_alpha10 = pnest),
  func.sor_pam10[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_alpha10= phylo.beta.SOR, func_turn_alpha10 = phylo.beta.SIM, func_nest_alpha10 = phylo.beta.SNE, func_pnest_alpha10 = pnest))


######################
## mapping spatial patterns of beta-diversity estimated ranges using different alpha hulls

# the template of raster for mapping 
raster_cells <- tree_pam6$Richness_Raster
raster_cells[] <- NA

# remove small islands in the shapefile of global lands
land_bhm <- land_bhm[land_bhm$area>500, ]

# the names of beta variables were wrote in the title
titles <- c("A Taxonomic total beta-diversity", "B Phylogenetic total beta-diversity", "C Functional total beta-diversity",
            "D Taxonomic turnover", "E Phylogenetic turnover", "F Functional turnover",
            "G Taxonomic nestedness", "H Phylogenetic nestedness", "I Functional nestedness",
            "J Taxonomic nestedness proportion", "K Phylogenetic nestedness proportion", "L Functional nestedness proportion")


###### beta diversity calculated with alpha hull 2
betaobs_alpha2 <- beta_alpha2[ , c(1,2,6,10,3,7,11,4,8,12,5,9,13)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(betaobs_alpha2$spe_total_alpha2, betaobs_alpha2$phylo_total_alpha2, betaobs_alpha2$func_total_alpha2)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(betaobs_alpha2$spe_turn_alpha2, betaobs_alpha2$phylo_turn_alpha2, betaobs_alpha2$func_turn_alpha2)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(betaobs_alpha2$spe_nest_alpha2, betaobs_alpha2$phylo_nest_alpha2, betaobs_alpha2$func_nest_alpha2)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(betaobs_alpha2$spe_pnest_alpha2, betaobs_alpha2$phylo_pnest_alpha2, betaobs_alpha2$func_pnest_alpha2)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)


# generate maps
beta.maps <- list()
for(i in c(1:12)){
  if(i %in% 1:3)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha2[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 4:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha2[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:9)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha2[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 10:12)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha2[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=3, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_Observed_beta_alpha2_maps.png", unit="mm", width=180, height=100)



###### beta diversity calculated with alpha hull 4
betaobs_alpha4 <- beta_alpha4[ , c(1,2,6,10,3,7,11,4,8,12,5,9,13)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(betaobs_alpha4$spe_total_alpha4, betaobs_alpha4$phylo_total_alpha4, betaobs_alpha4$func_total_alpha4)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(betaobs_alpha4$spe_turn_alpha4, betaobs_alpha4$phylo_turn_alpha4, betaobs_alpha4$func_turn_alpha4)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(betaobs_alpha4$spe_nest_alpha4, betaobs_alpha4$phylo_nest_alpha4, betaobs_alpha4$func_nest_alpha4)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(betaobs_alpha4$spe_pnest_alpha4, betaobs_alpha4$phylo_pnest_alpha4, betaobs_alpha4$func_pnest_alpha4)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)


# generate maps
beta.maps <- list()
for(i in c(1:12)){
  if(i %in% 1:3)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha4[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 4:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha4[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:9)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha4[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 10:12)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha4[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=3, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_Observed_beta_alpha4_maps.png", unit="mm", width=180, height=100)



###### beta diversity calculated with alpha hull 10
betaobs_alpha10 <- beta_alpha10[ , c(1,2,6,10,3,7,11,4,8,12,5,9,13)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(betaobs_alpha10$spe_total_alpha10, betaobs_alpha10$phylo_total_alpha10, betaobs_alpha10$func_total_alpha10)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(betaobs_alpha10$spe_turn_alpha10, betaobs_alpha10$phylo_turn_alpha10, betaobs_alpha10$func_turn_alpha10)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(betaobs_alpha10$spe_nest_alpha10, betaobs_alpha10$phylo_nest_alpha10, betaobs_alpha10$func_nest_alpha10)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(betaobs_alpha10$spe_pnest_alpha10, betaobs_alpha10$phylo_pnest_alpha10, betaobs_alpha10$func_pnest_alpha10)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)


# generate maps
beta.maps <- list()
for(i in c(1:12)){
  if(i %in% 1:3)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha10[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 4:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha10[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:9)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha10[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 10:12)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_alpha10[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=3, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_Observed_beta_alpha10_maps.png", unit="mm", width=180, height=100)



###################
# figure showing correlations between beta variablbes calculated with alpha hull 6 and other alpha hulls using different values of alphas

# set data.frames to store results
components <- c("Nestedness proportion","Nestedness","Turnover","Total")
dimensions <- c("Taxonommic","Phylogenetic","Functional")
beta.cor.ahull6.2 <- expand.grid(dimensions,components)
beta.cor.ahull6.4 <- expand.grid(dimensions,components)
beta.cor.ahull6.10 <- expand.grid(dimensions,components)
id <- c(5,9,13, 4,8,12, 3,7,11, 2,6,10)

# pearson's correlation between ahull6 with ahull2
match.rows <- intersect(beta_alpha6$ID, beta_alpha2$ID)
rows1 <- match(match.rows, beta_alpha6$ID)
rows2 <- match(match.rows, beta_alpha2$ID)
for(i in 1:12){
  beta.cor.ahull6.2[i,"cor"] <- round(cor(beta_alpha6[rows1, id[i]], beta_alpha2[rows2, id[i]]), 2)
}

# pearson's correlation between ahull6 with ahull4
match.rows <- intersect(beta_alpha6$ID, beta_alpha4$ID)
rows1 <- match(match.rows, beta_alpha6$ID)
rows2 <- match(match.rows, beta_alpha4$ID)
for(i in 1:12){
  beta.cor.ahull6.4[i,"cor"] <- round(cor(beta_alpha6[rows1, id[i]], beta_alpha4[rows2, id[i]]), 2)
}

# pearson's correlation between ahull6 with ahull10
match.rows <- intersect(beta_alpha6$ID, beta_alpha10$ID)
rows1 <- match(match.rows, beta_alpha6$ID)
rows2 <- match(match.rows, beta_alpha10$ID)
for(i in 1:12){
  beta.cor.ahull6.10[i,"cor"] <- round(cor(beta_alpha6[rows1, id[i]], beta_alpha10[rows2, id[i]]), 2)
}

# combine different comparison
beta.cor.ahull6.x <- bind_rows(beta.cor.ahull6.2 %>% mutate(panel = "Alpha 2 versus Alpha 6"),
          beta.cor.ahull6.4 %>% mutate(panel = "Alpha 4 versus Alpha 6"),
          beta.cor.ahull6.10 %>% mutate(panel = "Alpha 10 versus Alpha 6")) %>%
  mutate(panel = factor(panel, levels = c("Alpha 2 versus Alpha 6", "Alpha 4 versus Alpha 6", "Alpha 10 versus Alpha 6")))

# generate figure
library(egg) # a package used to generate tags
plot_beta_cor <- ggplot(data=beta.cor.ahull6.x, aes(Var1, Var2, fill= cor)) +
  facet_wrap(~panel) + 
  geom_tile(color="black", fill = "white") + 
  geom_text(aes(Var1, Var2, label = cor), color ="black", size = 2.8) +  
  labs(x = NULL, y = NULL) + 
  scale_y_discrete(labels = scales::wrap_format(4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=8, color="black"),
        axis.text.y = element_text(size=8, color="black", hjust = 0.5),
        panel.grid.major = element_blank())

# add tags
tag_facet(plot_beta_cor, vjust = 0, size = 2.81, open = "", close = "", tag_pool = LETTERS) + 
  coord_cartesian(clip = "off") + 
  theme(strip.text = element_text(size = 8))

ggsave(file="results/Fig.S_correlation_beta_betweenAlphas.png", unit="mm",width=180, height=65)

