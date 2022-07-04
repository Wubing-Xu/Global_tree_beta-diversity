###################################################################################################
# draw other supplementary figures and prepare supplementary tables

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
needed_libs <- c("tidyverse","ggplot2", "ggridges", "cowplot", "SpatialPack")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


source("Rcode/functions/scatter.plot.R")
source("Rcode/functions/mytmap.R")

load("models/beta_beta.deviation_envs.RDATA")
load("models/MM_SAR_beta_beta.deviation_envs.RDATA")
load("models/MM_SAR_beta_beta.deviation_envs_best.dist.RDATA")
load("intermediate_results/environments_final.RDATA")
load("data/tree_pam/tree_pam6_final.RDATA")
status_range <- read_csv("data/tree_pam/Status_Species_Ranges.csv")


####################
# scatter plot showing relationship between observed beta-diversity measures and temperature anomaly

# extract and organize needed variables 
beta_env_draw <- beta_env[, c(1, 5:8, 21, 9,13,17, 10,14,18, 11,15,19, 12,16,20)] %>% 
  pivot_longer(col =spe.beta.total:func.pnest, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("beta.", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  mutate(beta_facet = case_when(beta_facet == "spe" ~ "Taxonomic",
                                beta_facet == "phylo" ~ "Phylogenetic",
                                beta_facet == "func" ~ "Functional"),
         beta_facet = factor(beta_facet, levels = c("Taxonomic","Phylogenetic", "Functional")))

# prepare text added into figures
beta_anm_r2_text <- bind_rows(mm_sar_spe.total[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Taxonomic", beta_comp = "total"),
                              mm_sar_phylo.total[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Phylogenetic", beta_comp = "total"),
                              mm_sar_func.total[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Functional", beta_comp = "total"),
                              
                              mm_sar_spe.turn[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Taxonomic", beta_comp = "turn"),
                              mm_sar_phylo.turn[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Phylogenetic", beta_comp = "turn"),
                              mm_sar_func.turn[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Functional", beta_comp = "turn"),
                              
                              mm_sar_spe.nest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Taxonomic", beta_comp = "nest"),
                              mm_sar_phylo.nest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Phylogenetic", beta_comp = "nest"),
                              mm_sar_func.nest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Functional", beta_comp = "nest"),
                              
                              mm_sar_spe.pnest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Taxonomic", beta_comp = "pnest"),
                              mm_sar_phylo.pnest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Phylogenetic", beta_comp = "pnest"),
                              mm_sar_func.pnest[[1]][1, c("variable", "ps.cor", "cor.p")] %>% 
                                mutate(beta_facet = "Functional", beta_comp = "pnest"),
                              ) %>% 
  as_tibble() %>%
  mutate(beta_facet = factor(beta_facet, levels = c("Taxonomic", "Phylogenetic", "Functional")),
         label = LETTERS[1:12],
         r2 = ps.cor^2,
         r2_text = ifelse(cor.p<= 0.001, paste0('italic(R)^2 == ', round(r2, 3), '*"***"'), 
                          ifelse(cor.p<= 0.01 & cor.p> 0.001, paste0('italic(R)^2 == ', round(r2, 3), '*"**"'),
                                 ifelse(cor.p<= 0.05 & cor.p> 0.01, paste0('italic(R)^2 == ', round(r2, 3), '*"*"'),
                                        paste0('italic(R)^2 == ', round(r2, 3), "^ns")))))

# plot for total dissimilarity 
plot_total_anm <- ggplot(beta_env_draw %>% dplyr::filter(beta_comp == "total")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "total"), 
            aes(25, Inf, label = r2_text), hjust = 0.5, vjust = 2, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "total"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0.4, 0.99)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Total beta-diversity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for turnover 
plot_turn_anm <- ggplot(beta_env_draw %>% dplyr::filter(beta_comp == "turn")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "turn"), 
            aes(25, Inf, label = r2_text), hjust = 0.5, vjust = 2, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "turn"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0.1, 0.95)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Turnover") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for nestedness
plot_nest_anm <- ggplot(beta_env_draw %>% dplyr::filter(beta_comp == "nest")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "nest"), 
            aes(18, Inf, label = r2_text), hjust = 0.5, vjust = 2, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "nest"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0, 0.67)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Nestedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for nestedness proportion
plot_pnest_anm <- ggplot(beta_env_draw %>% dplyr::filter(beta_comp == "pnest")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(mat.anomaly, beta_values), alpha = 0.5) +
  geom_smooth(aes(mat.anomaly, beta_values), method = "lm", se = FALSE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "pnest"), 
            aes(18, Inf, label = r2_text), hjust = 0.5, vjust = 2, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_anm_r2_text  %>% dplyr::filter(beta_comp == "pnest"), 
            aes(3, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0, 0.94)) +
  labs(x = expression("Temperature anomaly ("*degree~C*")"), y = "Nestedness proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank())

cowplot::plot_grid(plot_total_anm, plot_turn_anm, plot_nest_anm, plot_pnest_anm,
                   nrow = 4, rel_heights =  c(1.15, 1, 1, 1.2))
ggsave(file="results/Fig.S_beta_anomaly.png", unit="mm",width=150, height=200)



##########################
## Scatter plot comparing species, functional and phylogenetic beta
beta_obs <- beta_env[, c(1, 5:8, 9,13,17, 10,14,18, 11,15,19, 12,16,20)] %>%
  as.data.frame()

# x- and y-labels
labels <- c("Taxonomic turnover", "Taxonomic nest", "Taxonomic total",
            "Phylogenetic turnover", "Phylogenetic nest", "Phylogenetic total",
            "Functional turnover", "Functional nest", "Functional total",
            "Taxonomic nest prop", "Phylogenetic nest prop", "Functional nest prop")

panel.labels <- LETTERS[1:12]

xid <- c(6,6,7, 9,9,10, 12,12,13, 15,15,16)
yid <- c(7,8,8, 10,10,11, 13,14,14, 16,17,17)
xlabs.id <- c(3,3,6,1,1,4,2,2,5,10,10,11)
ylabs.id <- c(6,9,9,4,7,7,5,8,8,11,12,12)

# generate figure
png("results/Fig.S_compare_beta_3dimentions.png", width=9, height=12, unit="in", res=300)
par(mfrow=c(4, 3), mar=c(5,5,1,1), cex=0.7)
for(i in 1:12){
  # determine extent of x and y 
  beta_xy <- c(beta_obs[,xid[i]], beta_obs[,yid[i]])
  beta_range <- max(beta_xy) - min(beta_xy)
  # plot
  scatter.plot(x=beta_obs[,xid[i]], y=beta_obs[,yid[i]], coords=beta_obs[,4:5], modified.ttest=TRUE, 
               xlim=c(min(beta_xy)-0.05*beta_range, max(beta_xy)+0.05*beta_range), 
               ylim=c(min(beta_xy)-0.05*beta_range, max(beta_xy)+0.05*beta_range),
               col=adjustcolor("black", alpha=0.2), col.line="red", 
               cex=0.6, cex.lab=1.5, cex.axis=1.3, text.pos=c(min(beta_xy)+0.35*beta_range, max(beta_xy)-0.1*beta_range), text.cex=1.5,
               xlab=labels[xlabs.id[i]], ylab=labels[ylabs.id[i]], pvalue.symbol=TRUE)
  mtext(side=2, panel.labels[i], line=3, las=1, at=max(beta_xy)+0.05*beta_range, cex=1.5)
  abline(0, 1, lty=2, col="blue", lwd=1.5)
}
dev.off()



####################
## get the importance values (weights), and correlations between environmental variables and observed beta-diversity measures,
# and show them as a figure and a table 

# extract weights and correlations from model resul tables
env_impor_beta.obs <- bind_rows(mm_sar_spe.total[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Total beta-diversity"),
          mm_sar_spe.turn[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Turnover"),
          mm_sar_spe.nest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Nestedness"),
          mm_sar_spe.pnest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Taxonomic", beta_comp = "Nestedness proportion"),
          
          mm_sar_phylo.total[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Total beta-diversity"),
          mm_sar_phylo.turn[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Turnover"),
          mm_sar_phylo.nest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Nestedness"),
          mm_sar_phylo.pnest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Phylogenetic", beta_comp = "Nestedness proportion"),
          
          mm_sar_func.total[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Total beta-diversity"),
          mm_sar_func.turn[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Turnover"),
          mm_sar_func.nest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Nestedness"),
          mm_sar_func.pnest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
            mutate(beta_facet = "Functional", beta_comp = "Nestedness proportion")
          ) %>%
  as_tibble()

#  format correlations and change names of environmental variables
env_impor_beta.obs <- env_impor_beta.obs %>%
  mutate(signif = ifelse(cor.p < 0.001, "***", 
                         ifelse(cor.p < 0.01 & cor.p > 0.001, "**", 
                                ifelse(cor.p < 0.05 & cor.p > 0.01, "*", ""))),
         cor_signif = paste0(round(ps.cor, 2), signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temp_anomaly",
                              variable == "map.anomaly" ~ "Prec_anomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temp_seas.",
                              variable == "ps" ~ "Prec_seas.",
                              variable == "log_topo" ~ "Elev_range")) %>%
  mutate(beta_facet = factor(beta_facet, levels = c("Taxonomic", "Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Total beta-diversity", "Turnover", "Nestedness", "Nestedness proportion")),
         variable = factor(variable, levels = c("Temp_anomaly", "Prec_anomaly", "MAT" ,"MAP", "Temp_seas.", "Prec_seas.", "Elev_range"))) 

# generate figure
beta_comp_color <- c("Total beta-diversity" = "#1B9E77", "Turnover" = "#7570B3", "Nestedness" = "#D95F02","Nestedness proportion" = "#E6AB02")

ggplot(env_impor_beta.obs) +
  facet_grid(beta_comp  ~ beta_facet) +
  geom_col(aes(y = fct_rev(variable), x = weight, fill = beta_comp)) + 
  labs(y = NULL, x = "Importance value (weights)") + 
  scale_fill_manual(name = NULL, values = beta_comp_color) + 
  theme_bw() +
  theme(legend.position = "no",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8))
  
ggsave(file="results/Fig.S_importance_envs_observed_beta.png", unit="mm",width=150, height=150)


## organize correlations between beta-diversity measures and environmental variables, 
# and save it as a table
table_cor_betaobs <- env_impor_beta.obs %>% 
  dplyr::select(beta_facet, beta_comp, variable, cor_signif) %>%
  pivot_wider(names_from = variable, values_from = cor_signif)

write.csv(table_cor_betaobs, file = "results/Table.S_correlation_obs.beta_envs.csv", row.names = FALSE)



###################
## get the importance values (weights), and correlations between environmental variables and 
# deviation of phylogenetic and functional turnover and nestedness; and show them as a figure and a table 

# extract weights and correlations from model resul tables
env_impor_beta.dev <- bind_rows(mm_sar_phylo.dev.turn[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
                                  mutate(beta_facet = "Phylogenetic", beta_comp = "Deviation of turnover"),
                                mm_sar_phylo.dev.nest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
                                  mutate(beta_facet = "Phylogenetic", beta_comp = "Deviation of nestedness"),
                                mm_sar_func.dev.turn[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
                                  mutate(beta_facet = "Functional", beta_comp = "Deviation of turnover"),
                                mm_sar_func.dev.nest[[1]][1:7, c("variable", "weight", "ps.cor", "cor.p")] %>% 
                                  mutate(beta_facet = "Functional", beta_comp = "Deviation of nestedness")
                                ) %>%
  as_tibble()

#  format correlations and change names of environmental variables
env_impor_beta.dev <- env_impor_beta.dev %>%
  mutate(signif = ifelse(cor.p < 0.001, "***", 
                         ifelse(cor.p < 0.01 & cor.p > 0.001, "**", 
                                ifelse(cor.p < 0.05 & cor.p > 0.01, "*", ""))),
         cor_signif = paste0(round(ps.cor, 2), signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temp_anomaly",
                              variable == "map.anomaly" ~ "Prec_anomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temp_seas.",
                              variable == "ps" ~ "Prec_seas.",
                              variable == "log_topo" ~ "Elev_range")) %>%
  mutate(beta_facet = factor(beta_facet, levels = c("Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Deviation of turnover", "Deviation of nestedness")),
         variable = factor(variable, levels = c("Temp_anomaly", "Prec_anomaly", "MAT" ,"MAP", "Temp_seas.", "Prec_seas.", "Elev_range"))) 

# generate figure
beta_comp_color <- c("Deviation of turnover" = "#7570B3", "Deviation of nestedness" = "#D95F02")

ggplot(env_impor_beta.dev) +
  facet_grid(beta_comp  ~ beta_facet) +
  geom_col(aes(y = fct_rev(variable), x = weight, fill = beta_comp)) + 
  labs(y = NULL, x = "Importance value (weights)") + 
  scale_fill_manual(name = NULL, values = beta_comp_color) + 
  theme_bw() +
  theme(legend.position = "no",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 8))

ggsave(file="results/Fig.S_importance_envs_beta_deviation.png", unit="mm",width=100, height=100)


## organize correlations between beta-diversity measures and environmental variables, 
# and save it as a table
table_cor_betadev <- env_impor_beta.dev %>% 
  dplyr::select(beta_facet, beta_comp, variable, cor_signif) %>%
  pivot_wider(names_from = variable, values_from = cor_signif)

write.csv(table_cor_betadev, file = "results/Table.S_correlation_beta.dev_envs.csv", row.names = FALSE)



#########################
## prepare the table showing best distance used in SAR-error models and related statistical variables

best_SAR_distance <- bind_rows(distSAR_spe.total[3, c(1, 5, 3, 6)],
          distSAR_spe.turn[3, c(1, 5, 3, 6)],
          distSAR_spe.nest[3, c(1, 5, 3, 6)],
          distSAR_spe.pnest[3, c(1, 5, 3, 6)],
          distSAR_phylo.total[3, c(1, 5, 3, 6)],
          distSAR_phylo.turn[3, c(1, 5, 3, 6)],
          distSAR_phylo.nest[3, c(1, 5, 3, 6)],
          distSAR_phylo.pnest[3, c(1, 5, 3, 6)],
          distSAR_func.total[3, c(1, 5, 3, 6)],
          distSAR_func.turn[3, c(1, 5, 3, 6)],
          distSAR_func.nest[3, c(1, 5, 3, 6)],
          distSAR_func.pnest[3, c(1, 5, 3, 6)],
          distSAR_phylo.dev.turn[3, c(1, 5, 3, 6)],
          distSAR_phylo.dev.nest[3, c(1, 5, 3, 6)],
          distSAR_func.dev.turn[3, c(1, 5, 3, 6)],
          distSAR_func.dev.nest[3, c(1, 5, 3, 6)]) %>%
  mutate(beta_facet = c(rep("Taxonomic", 4), rep("Phylogenetic", 4), rep("Functional", 4), rep("Phylogenetic", 2), rep("Functional", 2)),
         beta_comp = c(rep(c("Total dissimilarity", "Turnover", "Nestedness", "Nestedness proportion"), times = 3), 
                       rep(c("Deviation of turnover", "Deviation of nestedness"), times = 2))) %>% 
  mutate(across(c(I20, Imax, R2), round, 2)) %>%
  relocate(beta_facet, beta_comp)

write.csv(best_SAR_distance, file = "results/Table.S_best_SAR_distance.csv", row.names = FALSE)



####################
# maps for all environmental variables

# the template of raster for mapping 
raster_cells <- tree_pam6$Richness_Raster
raster_cells[] <- NA

# remove small islands in the shapefile of global lands
land_bhm <- land_bhm[land_bhm$area>500, ]

# the variables to map 
env200_2map <- beta_env %>%
  dplyr::select(ID, mat.anomaly, map.anomaly, mat, map, ts, ps, topo) %>%
  as.data.frame()

main.titles <- c("A Temperature anomaly since LGM", "B Precipitation anomaly since LGM", 
                 "C Mean annual temperature","D Mean annual precipitation",
                 "E Temperature seasonality","F Precipitation seasonality",
                 "G Elevation range")

# self_defined breaks for two map and elev_range
breaks.map <- c(0,800,1600,2400,3300)
breaks.eler <- c(0,1300,2600,3900,5200)

# generate maps
env.maps <- list()
for(i in c(1:7)){
  breaks <- NULL
  if(i==4) breaks <- breaks.map
  if(i==7) breaks <- breaks.eler
  env.map <- mytmap(myraster=raster_cells, myvalue=env200_2map[, c(1, i+1)], polygons=land_bhm, breaks=breaks, break.type="linear", mymidpoint="auto", 
                    digits=0, limit=0, cols=NULL, style="cont", 
                    main.title=main.titles[i], main.title.size=0.7, main.title.position=0.1,
                    legend.position=c(0.08, 0.01), legend.height=-0.6, legend.text.size=0.6)
  env.maps[[i]] <- env.map
}

# save figure
beta.maps <- tmap_arrange(env.maps, nrow=4, ncol=2, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_environmens_maps.png", unit="mm", width=120, height=100)


## map temperature anomaly since the LGM using all all grid cells (not just those with tree distributions) for preparing conceptual figure
mat_anomaly <- env200_mean_nn24_all[, c(1, 25)] %>% as.data.frame()
breaks.amomaly <- c(0, 7, 14, 21, 28, 35)
map_anomaly <- mytmap(myraster=raster_cells, myvalue=mat_anomaly[, c(1, 2)], polygons=land_bhm, breaks=breaks.amomaly, mymidpoint="auto", 
       digits=0, limit=0, cols=divPalette(100, "RdYlBu")[100:1], style="cont", 
       legend.position=c(0.08, 0.01), legend.height=-0.6, legend.text.size=0.6)
tmap_save(map_anomaly, file="results/Fig.S_temp_amomaly_map1.png", unit="mm", width=120, height=50, dpi = 1200)



####################
# the map for species richness

# calculate species richness
tree_sprich <- data.frame(ID = as.numeric(rownames(tree_pam6[[1]])), sprich = apply(tree_pam6[[1]][, -c(1:2)], 1, sum))

# generate the map
map_richness <- mytmap(myraster=raster_cells, myvalue=tree_sprich, polygons=land_bhm, breaks=NULL, break.type="log.linear", mymidpoint="auto", 
       digits=0, limit=0, cols=NULL, style="cont", 
       legend.position=c(0.08, 0.01), legend.height=-0.6, legend.text.size=0.6)

tmap_save(map_richness, file="results/Fig.S_species_richness.png", unit="mm", width=120, height=50, dpi = 300)



###########################
## figure showing correlation between environmental variables

# the used environmental variables
env200_2cor <- beta_env %>%
  dplyr::select(longitude, latitude, mat.anomaly, map.anomaly, mat, map, ts, ps, log_topo) %>%
  set_names(c("longitude", "latitude","Temp_anomaly", "Prec_anomaly", "MAT","MAP","Temp_seas.","Prec_seas.","Elev_range")) %>%
  as.data.frame()

# get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat,diag=TRUE)]<- NA
  return(cormat)
}

# get correlation betweenvariable pairs 
env200_cor <- round(cor(env200_2cor[, -(1:2)]),2) %>%
  get_upper_tri() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variable1") %>%
  pivot_longer(-1, names_to = "variable2", values_to = "cor_values") %>% 
  dplyr::filter(!variable1 == variable2 & !is.na(cor_values)) %>%
  mutate(pvalue = NA)

# add significance level using modified t-test
for(i in 1:nrow(env200_cor)){
  id1 <- which(colnames(env200_2cor) == env200_cor$variable1[i])
  id2 <- which(colnames(env200_2cor) == env200_cor$variable2[i])
  env200_cor$pvalue[i] <- modified.ttest(x=env200_2cor[, id1], y=env200_2cor[, id2], coords=env200_2cor[,1:2], nclass=NULL)$p.value
}

env200_cor <- env200_cor %>% 
  mutate(cor_label = ifelse(pvalue<= 0.001, paste0(cor_values, '*"***"'), 
                            ifelse(pvalue<= 0.01 & pvalue> 0.001, paste0(cor_values, '*"**"'),
                                   ifelse(pvalue<= 0.05 & pvalue> 0.01, paste0(cor_values, '*"*"'),
                                          cor_values)))) %>%
  mutate(variable1 = factor(variable1, levels = c("Temp_anomaly", "Prec_anomaly", "MAT","MAP","Temp_seas.","Prec_seas.","Elev_range")),
         variable2 = factor(variable2, levels = c("Temp_anomaly", "Prec_anomaly", "MAT","MAP","Temp_seas.","Prec_seas.","Elev_range")))

ggplot(data = env200_cor, aes(variable2, variable1,  fill = cor_values))+
  geom_tile(color = "black") +
  geom_text(aes(variable2, variable1, label = cor_label), color = "black", size = 2.81, parse =TRUE) + 
  labs(x = NULL, y = NULL) + 
  coord_fixed() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = c(0.4, 0.8),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 5, barheight = 0.8,
                               title.position = "top", title.hjust = 0.5))	

ggsave(file="results/Fig.S_correlation_environments.png", unit="mm", width=100, height=90)



###################
## figure showing latitudinal patterns of beta-diversity 
# fitted with segmented regression

# a self-defined function to fit segmented regression and get predicted values
library(segmented)
segmented_pred <- function(x, y){
  lm1 <- lm(y ~ x)
  seg1 <- segmented(lm1, ~ x, npsi = 1)
  y_pred <- predict(seg1)
  return(y_pred)
}

# extract and organize needed variables 
beta_lat <- beta_env[, c(1, 8, 9,13,17, 10,14,18, 11,15,19, 12,16,20)] %>% 
  pivot_longer(col =spe.beta.total:func.pnest, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("beta.", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  mutate(beta_facet = case_when(beta_facet == "spe" ~ "Taxonomic",
                                beta_facet == "phylo" ~ "Phylogenetic",
                                beta_facet == "func" ~ "Functional"),
         beta_facet = factor(beta_facet, levels = c("Taxonomic","Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("total","turn", "nest", "pnest")),
         latitude_abs = abs(latitude))

# regressing beta variables against latitude using segmented regressions and keep the predicted values
beta_lat <- beta_lat %>% 
  group_by(beta_comp, beta_facet) %>%
  mutate(beta_values_pred = segmented_pred(x = latitude_abs, y= beta_values)) 

beta_lat_r2 <- beta_lat %>%
  summarise(r2 = cor(beta_values_pred, beta_values)^2) %>% 
  ungroup() %>%
  mutate(r2_label = paste0('italic(R)^2 == ', round(r2, 3)),
         label = LETTERS[1:12])

# plot for total dissimilarity 
plot_total_lat <- ggplot(beta_lat %>% dplyr::filter(beta_comp == "total")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(latitude_abs, beta_values), alpha = 0.3) +
  geom_line(aes(latitude_abs, beta_values_pred), col = "blue", size = 1) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "total"), 
            aes(58, Inf, label = r2_label), hjust = 0.5, vjust = 1.5, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "total"), 
            aes(5, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0.4, 0.99)) +
  labs(x = expression("Absolute latitude ("*degree*")"), y = "Total beta-diversity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for turnover 
plot_turn_lat <- ggplot(beta_lat %>% dplyr::filter(beta_comp == "turn")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(latitude_abs, beta_values), alpha = 0.3) +
  geom_line(aes(latitude_abs, beta_values_pred), col = "blue", size = 1) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "turn"), 
            aes(58, Inf, label = r2_label), hjust = 0.5, vjust = 1.5, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "turn"), 
            aes(5, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0.1, 0.95)) +
  labs(x = expression("Absolute latitude ("*degree*")"), y = "Turnover") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for nestedness
plot_nest_lat <- ggplot(beta_lat %>% dplyr::filter(beta_comp == "nest")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(latitude_abs, beta_values), alpha = 0.3) +
  geom_line(aes(latitude_abs, beta_values_pred), col = "blue", size = 1) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "nest"), 
            aes(58, Inf, label = r2_label), hjust = 0.5, vjust = 1.5, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "nest"), 
            aes(5, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0, 0.67)) +
  labs(x = expression("Absolute latitude ("*degree*")"), y = "Nestedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

# plot for nestedness proportion
plot_pnest_lat <- ggplot(beta_lat %>% dplyr::filter(beta_comp == "pnest")) +
  facet_wrap( ~ beta_facet, scales = "fixed") +
  geom_point(aes(latitude_abs, beta_values), alpha = 0.3) +
  geom_line(aes(latitude_abs, beta_values_pred), col = "blue", size = 1) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "pnest"), 
            aes(58, Inf, label = r2_label), hjust = 0.5, vjust = 1.5, size = 2.81, parse =TRUE) + 
  geom_text(data = beta_lat_r2  %>% dplyr::filter(beta_comp == "pnest"), 
            aes(5, Inf, label = label), hjust = 0.5, vjust = 1.5, size = 3.515, fontface = "bold") + 
  scale_y_continuous(limits = c(0, 0.94)) +
  labs(x = expression("Absolute latitude ("*degree*")"), y = "Nestedness proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text  = element_blank())

cowplot::plot_grid(plot_total_lat, plot_turn_lat, plot_nest_lat, plot_pnest_lat,
                   nrow = 4, rel_heights =  c(1.15, 1, 1, 1.2))
ggsave(file="results/Fig.S_beta_latitude.png", unit="mm",width=150, height=200)



#############################
## figure showing histogram of the number of occurrences for each species 

# only keep species that were used in beta analyses
status_range1 <- status_range %>%
  dplyr::filter(sp %in% tree_pam6[[3]])

ggplot(status_range1) +
  geom_histogram(aes(N)) + 
  scale_x_log10() +
  labs(x = "Number of occurrences", y = "Number of species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"))

ggsave(file="results/Fig.S_histogram_nocc.png", unit="mm",width=80, height=80)
