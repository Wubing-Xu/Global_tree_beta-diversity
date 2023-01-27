###################################################################################################
## draw main figurs of this paper: fig.2 - fig. 6

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "XUWUBING-PC" = "D:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)


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
source("Rcode/functions/scatter.plot.R") 

load("models/beta_beta.deviation_envs.RDATA")
load("models/MM_SAR_beta_beta.deviation_envs.RDATA")
load("data/tree_pam/tree_pam6_final.RDATA")


####################
# figure 2: maps showing patterns of beta-diversity

# the template of raster for mapping 
raster_cells <- tree_pam6$Richness_Raster
raster_cells[] <- NA

# remove small islands in the shapefile of global lands
land_bhm <- land_bhm[land_bhm$area>500, ]

# the values to map 
beta_obs <- beta_env[, c(1,9,13,17 ,10,14,18,11,15,19,12,16,20)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(beta_obs$spe.beta.total, beta_obs$phylo.beta.total, beta_obs$func.beta.total)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(beta_obs$spe.beta.turn, beta_obs$phylo.beta.turn, beta_obs$func.beta.turn)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(beta_obs$spe.beta.nest, beta_obs$phylo.beta.nest, beta_obs$func.beta.nest)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(beta_obs$spe.pnest, beta_obs$phylo.pnest, beta_obs$func.pnest)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)

# the names of beta variables were wrote in the title
# titles <- c(expression("A TBD"[total]),expression("B PBD"[total]),expression("C FBD"[total]),
#             expression("D TBD"[turnover]),expression("E PBD"[turnover]),expression("F FBD"[turnover]),
#             expression("G TBD"[nestedness]),expression("H PBD"[nestedness]),expression("I FBD"[nestedness]),
#             expression("J TBD"[proportion.nest]),expression("K PBD"[proportion.nest]),expression("L FBD"[proportion.nest]))
titles <- c("A Taxonomic total beta-diversity", "B Phylogenetic total beta-diversity", "C Functional total beta-diversity",
            "D Taxonomic turnover", "E Phylogenetic turnover", "F Functional turnover",
            "G Taxonomic nestedness", "H Phylogenetic nestedness", "I Functional nestedness",
            "J Taxonomic nestedness proportion", "K Phylogenetic nestedness proportion", "L Functional nestedness proportion")

# generate maps
beta.maps <- list()
for(i in c(1:12)){
  if(i %in% 1:3)
    beta.map <- mytmap(myraster=raster_cells, myvalue=beta_obs[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.75, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 4:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=beta_obs[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.75, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:9)
    beta.map <- mytmap(myraster=raster_cells, myvalue=beta_obs[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.75, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 10:12)
    beta.map <- mytmap(myraster=raster_cells, myvalue=beta_obs[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.75, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=3, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.2_Observed_beta_maps.png", unit="mm", width=180, height=100)
tmap_save(beta.maps, file="results/Fig.2_Observed_beta_maps.pdf", unit="mm", width=180, height=100)



####################
# figure 4: maps showing patterns of residuals of observed phylogenetic and functional turnover and nestedness from the random expected
# the values to map 
beta_dev <- betadev_env[, c(1, 10, 14, 11, 15)] %>%
  as.data.frame()

# the names of variable were wrote in the title
#titles <- c(expression("A Deviation of PBD"[turnover]), expression("B Deviation of FBD"[turnover]),
#             expression("C Deviation of PBD"[nestedness]), expression("D Deviation of FBD"[nestedness]))
titles <- c("A Deviation of phylogenetic turnover", "B Deviation of functional turnover",
            "C Deviation of phylogenetic nestedness", "D Deviation of functional nestedness")


# generate maps
beta.dev.maps <- list()
for(i in c(1:4)){
  beta.map <- mytmap(myraster=raster_cells, myvalue=beta_dev[,c(1, i+1)], polygons=land_bhm, breaks=NULL, break.type="linear", mymidpoint=0, 
                     digits=2, limit=0.01, cols=NULL, style="cont", inner.margins=c(0.01,0,0.15,0), 
                     title=titles[i], title.size=0.75, title.position = c(0.10, 0.92),
                     legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.dev.maps[[i]] <- beta.map
}

beta.dev.maps <- tmap_arrange(beta.dev.maps, nrow=2, ncol=2, outer.margins=0)
tmap_save(beta.dev.maps, file="results/Fig.4_Deviation_PBD.FBD_maps.png", unit="mm",width=150, height=67)
tmap_save(beta.dev.maps, file="results/Fig.4_Deviation_PBD.FBD_maps.pdf", unit="mm",width=150, height=67)



####################
# figure 3: plots showing standardized coefficients from models regressing beta-diversity variables against environmental variables 

# extract standardized averaged coefficients, their 95% CI and p-values from models
mmsar_coef_beta.obs <- bind_rows(mm_sar_spe.total[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Total beta-diversity"),
          mm_sar_spe.turn[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Turnover"),
          mm_sar_spe.nest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Nestedness"),
          mm_sar_spe.pnest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Nestedness proportion"),
          
          mm_sar_phylo.total[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Total beta-diversity"),
          mm_sar_phylo.turn[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Turnover"),
          mm_sar_phylo.nest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Nestedness"),
          mm_sar_phylo.pnest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Nestedness proportion"),
          
          mm_sar_func.total[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Total beta-diversity"),
          mm_sar_func.turn[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Turnover"),
          mm_sar_func.nest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Nestedness"),
          mm_sar_func.pnest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Nestedness proportion")
          ) %>%
  as_tibble()

# add significance and changes names of environmental variables
mmsar_coef_beta.obs <- mmsar_coef_beta.obs %>%
  mutate(signif = ifelse(p.avg.cod < 0.05, "1", "0"),
         beta_signif = ifelse(signif==1, beta_comp, signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temperature\nanomaly",
                              variable == "map.anomaly" ~ "Precipitation\nanomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temperature\nseasonality",
                              variable == "ps" ~ "Precipitation\nseasonality",
                              variable == "log_topo" ~ "Elevational\nrange", 
                              variable == "sqrt_hmi" ~ "Human\nmodification"))

# set factor levels
mmsar_coef_beta.obs <- mmsar_coef_beta.obs %>% 
  mutate(beta_facet = factor(beta_facet, levels = c("Taxonomic", "Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Total beta-diversity", "Turnover", "Nestedness", "Nestedness proportion")),
         variable = factor(variable, levels = c("Temperature\nanomaly", "Precipitation\nanomaly", "MAT" ,"MAP", "Temperature\nseasonality", 
                                                "Precipitation\nseasonality", "Elevational\nrange", "Human\nmodification"))) 

# shape and colours for components of beta; colours of all components that are not significant are set as gray
beta_comp_shape <- c("Total beta-diversity" = 21, "Turnover" = 22, "Nestedness" = 23, "Nestedness proportion" = 24)
# beta_comp_color <- c("Total" = 2, "Turnover" = 8, "Nestedness" = 4, "Prop.nest" = 5, "0" = "gray")
beta_comp_color <- c("Total beta-diversity" = "#1B9E77", "Turnover" = "#7570B3", "Nestedness" = "#D95F02","Nestedness proportion" = "#E6AB02", "0" = "gray")

## generate figure
# get legend showing different components of beta
plot_legend_beta_comp <- ggplot(mmsar_coef_beta.obs) + 
  facet_wrap(~ beta_facet, scale = "free_x") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(beta_comp), color = beta_comp), 
                  position = position_dodge(0.6), size =0.8) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(beta_comp), fill = beta_comp, color = beta_comp),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  labs(y = NULL, x = "Standardized Coefficient") +  
  scale_shape_manual(name = NULL, values = beta_comp_shape) + 
  scale_color_manual(name = NULL, values = beta_comp_color) + 
  scale_fill_manual(name = NULL, values = beta_comp_color) + 
  guides(shape = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(legend.position = "bottom", # "bottom", c(0.92, 0.18)
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.text = element_text(size = 9))

legend_beta_comp <- get_legend(plot_legend_beta_comp)

# figure
plot_coef_obsbeta <- ggplot(mmsar_coef_beta.obs) + 
  facet_wrap(~ beta_facet, scale = "fixed") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(beta_comp), color = beta_signif), 
                  position = position_dodge(0.6), size =0.8) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(beta_comp), fill = beta_signif, color = beta_signif),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  scale_x_continuous(limit = c(-1.05, 1.05)) +
  labs(y = NULL, x = "Standardized coefficient") +  
  scale_shape_manual(name = NULL, values = beta_comp_shape) + 
  scale_color_manual(name = NULL, values = beta_comp_color, guide ="none") + 
  scale_fill_manual(name = NULL, values = beta_comp_color, guide ="none") + 
  theme_bw() +
  theme(legend.position = "n", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,0,1), "mm"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 9))

cowplot::plot_grid(plot_coef_obsbeta, plot_grid(NULL, legend_beta_comp, rel_widths =  c(0.15, 1)),
                   nrow = 2,  rel_heights =  c(1, 0.1))
ggsave(file="results/Fig.3_coefficients_observed_beta.png", unit="mm",width=180, height=150)
ggsave(file="results/Fig.3_coefficients_observed_beta-1.pdf", unit="mm",width=180, height=150)


####################
# figure 5: plots showing standardized coefficients from models regressing deviation of phylogenetic and functional
# turnover and nestedness against environmental variables 

# extract standardized averaged coefficients, their 95% CI and p-values from models
mmsar_coef_beta.dev <- bind_rows(mm_sar_phylo.dev.turn[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
                                    mutate(beta_facet = "Phylogenetic", beta_comp = "Deviation of turnover"),
                                  mm_sar_phylo.dev.nest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
                                    mutate(beta_facet = "Phylogenetic", beta_comp = "Deviation of nestedness"),
                                  mm_sar_func.dev.turn[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
                                    mutate(beta_facet = "Functional", beta_comp = "Deviation of turnover"),
                                  mm_sar_func.dev.nest[[1]][1:8, c("variable", "coef.avg.cod", "coef_2.5", "coef_97.5", "p.avg.cod")] %>% 
                                    mutate(beta_facet = "Functional", beta_comp = "Deviation of nestedness")) %>% 
  as_tibble()

# add significance and changes names of environmental variables
mmsar_coef_beta.dev <- mmsar_coef_beta.dev %>%
  mutate(signif = ifelse(p.avg.cod < 0.05, "1", "0"),
         beta_signif = ifelse(signif==1, beta_comp, signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temperature\nanomaly",
                              variable == "map.anomaly" ~ "Precipitation\nanomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temperature\nseasonality",
                              variable == "ps" ~ "Precipitation\nseasonality",
                              variable == "log_topo" ~ "Elevational\nrange", 
                              variable == "sqrt_hmi" ~ "Human\nmodification"))

# set factor levels
mmsar_coef_beta.dev <- mmsar_coef_beta.dev %>% 
  mutate(beta_facet = factor(beta_facet, levels = c("Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Deviation of turnover", "Deviation of nestedness")),
         variable = factor(variable, levels = c("Temperature\nanomaly", "Precipitation\nanomaly", "MAT" ,"MAP", "Temperature\nseasonality", 
                                                "Precipitation\nseasonality", "Elevational\nrange", "Human\nmodification"))) 

# shape and colours for components of beta; colours of all components that are not significant are set as gray
beta_comp_shape <- c("Deviation of turnover" = 22, "Deviation of nestedness" = 23)
# beta_comp_color <- c("Total" = 2, "Turnover" = 8, "Nestedness" = 4, "Prop.nest" = 5, "0" = "gray")
beta_comp_color <- c("Deviation of turnover" = "#7570B3", "Deviation of nestedness" = "#D95F02", "0" = "gray")

## generate figure
# get legend showing different components of beta
plot_legend_betadev_comp <- ggplot(mmsar_coef_beta.dev) + 
  facet_wrap(~ beta_facet, scale = "free_x") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(beta_comp), color = beta_comp), 
                  position = position_dodge(0.6), size =0.8) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(beta_comp), fill = beta_comp, color = beta_comp),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  labs(y = NULL, x = "Standardized Coefficient") +  
  scale_shape_manual(name = NULL, values = beta_comp_shape) + 
  scale_color_manual(name = NULL, values = beta_comp_color) + 
  scale_fill_manual(name = NULL, values = beta_comp_color) + 
  guides(shape = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(legend.position = "bottom", # "bottom", c(0.92, 0.18)
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 0, 0))

legend_betadev_comp <- get_legend(plot_legend_betadev_comp)

# figure
plot_coef_betadev <- ggplot(mmsar_coef_beta.dev) + 
  facet_wrap(~ beta_facet, scale = "fixed") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(beta_comp), color = beta_signif), 
                  position = position_dodge(0.6), size =0.8) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(beta_comp), fill = beta_signif, color = beta_signif),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  scale_x_continuous(limit = c(-1.35, 1.35)) +
  labs(y = NULL, x = "Standardized coefficient") +  
  scale_shape_manual(name = NULL, values = beta_comp_shape) + 
  scale_color_manual(name = NULL, values = beta_comp_color, guide ="none") + 
  scale_fill_manual(name = NULL, values = beta_comp_color, guide ="none") + 
  theme_bw() +
  theme(legend.position = "n", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,0,1), "mm"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 9),
        legend.text = element_text(size = 9))

cowplot::plot_grid(plot_coef_betadev, plot_grid(NULL, legend_betadev_comp, rel_widths =  c(0.2, 1)),
                   nrow = 2,  rel_heights =  c(1, 0.1))
ggsave(file="results/Fig.5_coefficients_beta_deviation.png", unit="mm",width=127, height=120)
ggsave(file="results/Fig.5_coefficients_beta_deviation.pdf", unit="mm",width=127, height=120)



####################
# figure 6: scatter plot showing relationship between deviation of turnover and nestedness and temperature anomaly

# extract and organive needed varialbes 
betadev_env_draw <- betadev_env[, c(1, 5:8, 17, 10, 11, 14, 15)] %>% 
  #rename(Phylogenetic_Turnover = phylo.beta.turn, phylo.beta.nest, func.beta.turn, func.beta.nest)
  pivot_longer(col =phylo.beta.turn:func.beta.nest, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("beta.", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  mutate(beta_facet = case_when(beta_facet == "phylo" ~ "Phylogenetic",
                                beta_facet == "func" ~ "Functional"),
         beta_facet = factor(beta_facet, levels = c("Phylogenetic", "Functional")))

# prepare text added into figures
betadev_anm_r2_text <- bind_rows(mm_sar_phylo.dev.turn[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                   mutate(beta_facet = "Phylogenetic", beta_comp = "turn"),
                                 mm_sar_func.dev.turn[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                   mutate(beta_facet = "Functional", beta_comp = "turn"),
                                 mm_sar_phylo.dev.nest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                   mutate(beta_facet = "Phylogenetic", beta_comp = "nest"),
                                 mm_sar_func.dev.nest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                   mutate(beta_facet = "Functional", beta_comp = "nest")) %>% 
  as_tibble() %>%
  mutate(beta_facet = factor(beta_facet, levels = c("Phylogenetic", "Functional")),
         label = c("A", "B", "C", "D"),
         r2 = ps.cor^2,
         r2_text = ifelse(cor.p<= 0.001, paste0('italic(R)^2 == ', round(r2, 3), '*"***"'), 
                          ifelse(cor.p<= 0.01 & cor.p> 0.001, paste0('italic(R)^2 == ', round(r2, 3), '*"**"'),
                                 ifelse(cor.p<= 0.05 & cor.p> 0.01, paste0('italic(R)^2 == ', round(r2, 3), '*"*"'),
                                        paste0('italic(R)^2 == ', round(r2, 3), "^ns")))))

# plot for deviation of turnover 
plot_devturn_anm <- ggplot(betadev_env_draw %>% dplyr::filter(beta_comp == "turn")) +
  facet_wrap( ~ beta_facet, scales = "free_y") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = betadev_anm_r2_text  %>% dplyr::filter(beta_comp == "turn"), 
            aes(25, Inf, label = r2_text), hjust = 0.5, vjust = 1.5, size = 2.81, parse =TRUE) + 
  geom_text(data = betadev_anm_r2_text  %>% dplyr::filter(beta_comp == "turn"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  #scale_y_continuous(limits = c(-0.3, 0.15)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Deviation of turnover") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for deviation of nestedness
plot_devnest_anm <- ggplot(betadev_env_draw %>% dplyr::filter(beta_comp == "nest")) +
  facet_wrap( ~ beta_facet, scales = "free_y") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = betadev_anm_r2_text  %>% dplyr::filter(beta_comp == "nest"), 
            aes(18, Inf, label = r2_text), hjust = 0.5, vjust = 2, size = 2.81, parse =TRUE) + 
  geom_text(data = betadev_anm_r2_text  %>% dplyr::filter(beta_comp == "nest"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  #scale_y_continuous(limits = c(-0.1, 0.3)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Devation of nestedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank(),)

cowplot::plot_grid(plot_devturn_anm, plot_devnest_anm, 
                   nrow = 2, rel_heights =  c(1, 1.05))
ggsave(file="results/Fig.6_beta_deviation_anomaly.png", unit="mm",width=120, height=110)
ggsave(file="results/Fig.6_beta_deviation_anomaly.pdf", unit="mm",width=120, height=110)
