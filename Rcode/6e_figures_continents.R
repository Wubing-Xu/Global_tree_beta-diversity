###################################################################################################
## draw figures showing standardized coefficients of environmental variables across continents
# and figures comparing phylogenetic and functional beta-diversity for each of continents

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


load("models/MM_SAR_region_beta_beta.deviation_envs.RDATA")


####################
## plots showing standardized coefficients from models regressing beta-diversity variables against environmental variables for each continent

# extract standardized averaged coefficients, their 95% CI and p-values from models
continents <- rep(names(mm_sar_region_spe.total), each = nrow(mm_sar_region_spe.total[[1]][[1]]))

mmsar_coef_region_beta.obs <- bind_rows(do.call(bind_rows, lapply(mm_sar_region_spe.total, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Taxonomic", beta_comp = "Total beta-diversity"),
          do.call(bind_rows, lapply(mm_sar_region_spe.turn, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Taxonomic", beta_comp = "Turnover"),
          do.call(bind_rows, lapply(mm_sar_region_spe.nest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Taxonomic", beta_comp = "Nestedness"),
          do.call(bind_rows, lapply(mm_sar_region_spe.pnest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Taxonomic", beta_comp = "Nestedness proportion"),
          
          do.call(bind_rows, lapply(mm_sar_region_phylo.total, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Total beta-diversity"),
          do.call(bind_rows, lapply(mm_sar_region_phylo.turn, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Turnover"),
          do.call(bind_rows, lapply(mm_sar_region_phylo.nest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Nestedness"),
          do.call(bind_rows, lapply(mm_sar_region_phylo.pnest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Nestedness proportion"),
          
          do.call(bind_rows, lapply(mm_sar_region_func.total, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Functional", beta_comp = "Total beta-diversity"),
          do.call(bind_rows, lapply(mm_sar_region_func.turn, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Functional", beta_comp = "Turnover"),
          do.call(bind_rows, lapply(mm_sar_region_func.nest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Functional", beta_comp = "Nestedness"),
          do.call(bind_rows, lapply(mm_sar_region_func.pnest, "[[", 1)) %>% 
            mutate(continent = continents, beta_facet = "Functional", beta_comp = "Nestedness proportion")) %>%
  as_tibble() %>%
  dplyr::select(variable, coef.avg.cod, coef_2.5, coef_97.5, p.avg.cod, continent, beta_facet, beta_comp) %>%
  dplyr::filter(!is.na(coef.avg.cod))




# add significance and changes names of environmental variables
mmsar_coef_region_beta.obs <- mmsar_coef_region_beta.obs %>%
  mutate(signif = ifelse(p.avg.cod < 0.05, "1", "0"),
         continent_signif = ifelse(signif==1, continent, signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temp_anomaly",
                              variable == "map.anomaly" ~ "Prec_anomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temp_seas.",
                              variable == "ps" ~ "Prec_seas.",
                              variable == "log_topo" ~ "Elev_range", 
                              variable == "sqrt_hmi" ~ "Human_modif."))

# set factor levels
mmsar_coef_region_beta.obs <- mmsar_coef_region_beta.obs %>% 
  mutate(beta_facet = factor(beta_facet, levels = c("Taxonomic", "Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Total beta-diversity", "Turnover", "Nestedness", "Nestedness proportion")),
         variable = factor(variable, levels = c("Temp_anomaly", "Prec_anomaly", "MAT" ,"MAP", "Temp_seas.", "Prec_seas.", "Elev_range", "Human_modif."))) 


library(RColorBrewer)
display.brewer.all(n=6, type="qual", select="Dark2", exact.n=TRUE)
brewer.pal.info(n=6, type="qual", select="Dark2", exact.n=TRUE)              
brewer.pal(6, "Dark2")

# shape and colours for continents; colours of coefficients that are not significant are set as gray
continent_shape <- c("Africa" = 21, "Asia" = 22 , "Australia" = 23 ,"Europe" = 24,"North America" = 4, "South America" = 8)
continent_color <- c("Africa" = "#1B9E77", "Asia" = "#D95F02" , "Australia" =  "#7570B3" ,"Europe" = "#E7298A", 
                     "North America" = "#66A61E", "South America" = "#E6AB02", "0" = "gray")

## generate figure
# get legend showing different continents
plot_legend_continent <- ggplot(mmsar_coef_region_beta.obs) + 
  facet_grid(beta_comp ~ beta_facet, scale = "fixed") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(continent), color = continent), 
                  position = position_dodge(0.6), size =0.4) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(continent), fill = continent, color = continent),
             position = position_dodge(0.6), size = 1) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  # scale_x_continuous(limit = c(-1.05, 1.05)) +
  labs(y = NULL, x = "Standardized coefficient") +  
  scale_shape_manual(name = NULL, values = continent_shape) + 
  scale_color_manual(name = NULL, values = continent_color) + 
  scale_fill_manual(name = NULL, values = continent_color) + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.margin = margin(0,10,0,10),
        legend.box.spacing = unit(0, units = 'mm'),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.text = element_text(size = 9)) +
  guides(shape = guide_legend(nrow = 1, reverse = TRUE)) +
  guides(color = guide_legend(nrow = 1)) +
  guides(fill = guide_legend(nrow = 1))

legend_continent <- get_legend(plot_legend_continent)

plot_coef_obsbeta <- ggplot(mmsar_coef_region_beta.obs) + 
  facet_grid(beta_comp ~ beta_facet, scale = "fixed") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(continent), color = continent_signif), 
                  position = position_dodge(0.6), size =0.4) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(continent), fill = continent_signif, color = continent_signif),
             position = position_dodge(0.6), size = 1) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  # scale_x_continuous(limit = c(-1.05, 1.05)) +
  labs(y = NULL, x = "Standardized coefficient") +  
  scale_shape_manual(name = NULL, values = continent_shape) + 
  scale_color_manual(name = NULL, values = continent_color, guide ="none") + 
  scale_fill_manual(name = NULL, values = continent_color, guide ="none") + 
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

cowplot::plot_grid(plot_coef_obsbeta, legend_continent,
                   nrow = 2, rel_heights =  c(1, 0.03))

ggsave(file="results/Fig.S_coefficients_observed_beta_continents.png", unit="mm",width=180, height=230)



####################
# plots showing standardized coefficients from models regressing deviation of phylogenetic and functional
# turnover and nestedness against environmental variables for each continent

# extract standardized averaged coefficients, their 95% CI and p-values from models
continents <- rep(names(mm_sar_region_phylo.dev.turn), each = nrow(mm_sar_region_phylo.dev.turn[[1]][[1]]))

mmsar_coef_region_beta.dev <- bind_rows(do.call(bind_rows, lapply(mm_sar_region_phylo.dev.turn, "[[", 1)) %>% 
                                          mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Deviation of turnover"),
                                        do.call(bind_rows, lapply(mm_sar_region_phylo.dev.nest, "[[", 1)) %>% 
                                          mutate(continent = continents, beta_facet = "Phylogenetic", beta_comp = "Deviation of nestedness"),
                                        do.call(bind_rows, lapply(mm_sar_region_func.dev.turn, "[[", 1)) %>% 
                                          mutate(continent = continents, beta_facet = "Functional", beta_comp = "Deviation of turnover"),
                                        do.call(bind_rows, lapply(mm_sar_region_func.dev.nest, "[[", 1)) %>% 
                                          mutate(continent = continents, beta_facet = "Functional", beta_comp = "Deviation of nestedness")) %>%
  as_tibble() %>%
  dplyr::select(variable, coef.avg.cod, coef_2.5, coef_97.5, p.avg.cod, continent, beta_facet, beta_comp) %>%
  dplyr::filter(!is.na(coef.avg.cod))


# add significance and changes names of environmental variables
mmsar_coef_region_beta.dev <- mmsar_coef_region_beta.dev %>%
  mutate(signif = ifelse(p.avg.cod < 0.05, "1", "0"),
         continent_signif = ifelse(signif==1, continent, signif),
         variable = case_when(variable == "mat.anomaly" ~ "Temp_anomaly",
                              variable == "map.anomaly" ~ "Prec_anomaly",
                              variable == "mat" ~ "MAT",
                              variable == "map" ~ "MAP",
                              variable == "ts" ~ "Temp_seas.",
                              variable == "ps" ~ "Prec_seas.",
                              variable == "log_topo" ~ "Elev_range", 
                              variable == "sqrt_hmi" ~ "Human_modif."))

# set factor levels
mmsar_coef_region_beta.dev <- mmsar_coef_region_beta.dev %>% 
  mutate(beta_facet = factor(beta_facet, levels = c("Phylogenetic", "Functional")),
         beta_comp = factor(beta_comp, levels = c("Deviation of turnover", "Deviation of nestedness")),
         variable = factor(variable, levels = c("Temp_anomaly", "Prec_anomaly", "MAT" ,"MAP", "Temp_seas.", "Prec_seas.", "Elev_range", "Human_modif."))) 


## generate figure
plot_coef_betadev <- ggplot(mmsar_coef_region_beta.dev) + 
  facet_grid(beta_comp ~ beta_facet, scale = "fixed") + 
  geom_linerangeh(aes(xmin=coef_2.5, xmax=coef_97.5, y = fct_rev(variable), group = fct_rev(continent), color = continent_signif), 
                  position = position_dodge(0.6), size =0.4) +
  geom_point(aes(y = variable, x = coef.avg.cod, shape = fct_rev(continent), fill = continent_signif, color = continent_signif),
             position = position_dodge(0.6), size = 1) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  # scale_x_continuous(limit = c(-1.05, 1.05)) +
  labs(y = NULL, x = "Standardized coefficient") +  
  scale_shape_manual(name = NULL, values = continent_shape) + 
  scale_color_manual(name = NULL, values = continent_color, guide ="none") + 
  scale_fill_manual(name = NULL, values = continent_color, guide ="none") + 
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

cowplot::plot_grid(plot_coef_betadev, legend_continent,
                   nrow = 2, rel_heights =  c(1, 0.03))

ggsave(file="results/Fig.S_coefficients_beta_deviation_continents.png", unit="mm",width=140, height=150)



####################################
## Scatter plot comparing functional and phylogenetic beta diversity for each of six continents

## organize phylogenetic and functional beta diversity variables
phylo_func_beta <- beta_env %>%
  dplyr::select(ID, continent, longitude, latitude, phylo.beta.total:func.pnest) %>%
  pivot_longer(col = 5:12, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("beta.", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  pivot_wider(names_from = beta_facet, values_from = beta_values) %>%
  mutate(beta_comp = case_when(beta_comp == "total" ~ "Total beta-diversity",
                               beta_comp == "turn" ~ "Turnover",
                               beta_comp == "nest" ~ "Nestedness",
                               beta_comp == "pnest" ~ "Nestedness proportion"),
         beta_comp = factor(beta_comp, levels = c("Total beta-diversity", "Turnover", "Nestedness", "Nestedness proportion")))

# calculate correlation between phylognetic and functional beta-diversity and significance using modified t-test
betacomp_continent <- phylo_func_beta %>% distinct(beta_comp, continent) %>%
  mutate(cor.r = NA,
         cor.p = NA)

for(i in 1:nrow(betacomp_continent)){
  phylo_func_beta_sub <- phylo_func_beta %>% dplyr::filter(beta_comp == betacomp_continent$beta_comp[i] & continent == betacomp_continent$continent[i])
  xy.cor <- modified.ttest(x=phylo_func_beta_sub$phylo, y=phylo_func_beta_sub$func, coords=phylo_func_beta_sub[, c("longitude", "latitude")], nclass=NULL)
  betacomp_continent$cor.r[i] = xy.cor$corr
  betacomp_continent$cor.p[i] = xy.cor$p.value
}

# prepare correlation text added into figures
betacomp_continent <- betacomp_continent %>%
  mutate(r_text = ifelse(cor.p<= 0.001, paste0('italic(r) == ', round(cor.r, 3), '*"***"'), 
                         ifelse(cor.p<= 0.01 & cor.p> 0.001, paste0('italic(r) == ', round(cor.r, 3), '*"**"'),
                                ifelse(cor.p<= 0.05 & cor.p> 0.01, paste0('italic(r) == ', round(cor.r, 3), '*"*"'),
                                       paste0('italic(r) == ', round(cor.r, 3), "^ns")))))

# generate scatter figures 
ggplot(phylo_func_beta) +
  facet_grid(continent ~ beta_comp, scales = "fixed") +
  geom_point(aes(phylo, func), alpha = 0.2, size = 0.4) +
  geom_smooth(aes(phylo, func), method = "lm", se = FALSE, color = "red", size = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = "blue", lty = 2, size = 0.5) +
  geom_text(data = betacomp_continent, 
            aes(0.4, Inf, label = r_text), hjust = 0.5, vjust = 2.5, size = 2.81, parse =TRUE) + 
  coord_fixed() +
  scale_x_continuous(limits = c(0, 0.95), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_y_continuous(limits = c(0, 0.95), breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(x = "Phylogenetic", y = "Functional") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 9))

ggsave(file="results/Fig.S_compare_phylo_func_beta_continents.png.png", unit="mm",width=180, height=260, dpi = 600)

