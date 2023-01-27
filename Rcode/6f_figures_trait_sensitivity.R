###################################################################################################
# prepare supplementary figures and tables for sensitivity analyses of trait imputation 

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "XUWUBING-PC" = "D:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)


# load packages
needed_libs <- c("tidyverse","ggplot2", "ggridges", "cowplot", "SpatialPack", "phytools", "parallel")

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

load("intermediate_results/observed_beta_sorensen.RDATA")
load("intermediate_results/observed_beta_sorensen_SpeciesWithTraits.RDATA")
load("models/beta_beta.deviation_envs.RDATA")
load("data/tree_pam/tree_pam6_final.RDATA")
load("data/traits_sensitivity/traits_coverage_sensitivity.RDATA")
load("data/traits_phylogeny/trait_phylogeny_final.RDATA")


##############################
## Map showing trait coverage (proportion of species with trait data)

# the template of raster for mapping 
raster_cells <- tree_pam6$Richness_Raster
raster_cells[] <- NA

# remove small islands in the shapefile of global lands
land_bhm <- land_bhm[land_bhm$area>500, ]

# the variables to map 
trait_coverage_2map <- trait_coverage %>%
  dplyr::select(ID, prop_spp_1trait, prop_spp_3trait, prop_spp_5trait) %>%
  as.data.frame()

trait_coverage_2map <- trait_coverage %>%
  as.data.frame()


trait1_coverage <- mytmap(myraster=raster_cells, myvalue=trait_coverage_2map[, c(1, 3)], polygons=land_bhm, breaks=NULL, break.type="linear", mymidpoint="auto", 
                          digits=2, limit=0, cols=NULL, style="cont", 
                          main.title="A Proportion of species with at least one trait", main.title.size=0.7, main.title.position = 0.1,
                          legend.position=c(0.07, 0.01), legend.height=-0.6, legend.text.size=0.6)

trait5_coverage <- mytmap(myraster=raster_cells, myvalue=trait_coverage_2map[, c(1, 7)], polygons=land_bhm, breaks=NULL, break.type="linear", mymidpoint="auto", 
                          digits=2, limit=0, cols=NULL, style="cont", 
                          main.title="B Proportion of species with at least five traits", main.title.size=0.7, main.title.position = 0.1,
                          legend.position=c(0.07, 0.01), legend.height=-0.6, legend.text.size=0.6)

trait8_coverage <- mytmap(myraster=raster_cells, myvalue=trait_coverage_2map[, c(1, 10)], polygons=land_bhm, breaks=NULL, break.type="linear", mymidpoint="auto", 
                          digits=2, limit=0, cols=NULL, style="cont", 
                          main.title="C Proportion of species with all eight traits", main.title.size=0.7, main.title.position = 0.1,
                          legend.position=c(0.07, 0.01), legend.height=-0.6, legend.text.size=0.6)

maps_trait_coverage <- tmap_arrange(trait1_coverage, trait5_coverage, trait8_coverage, 
                              nrow=3, ncol=1, outer.margins=0)
tmap_save(maps_trait_coverage, file="results/Fig.S_trait_coverage1.png",unit="mm",width=90, height=120)



##############################
## mapping spatial patterns of beta-diversity using species with at least one or five traits


## organize beta diversity variable using species with at least one or five traits, or all species
# at least one trait
beta_1trait <- bind_cols(
  phylo.sor_pam1trait[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_1trait= phylo.beta.SOR, phylo_turn_1trait = phylo.beta.SIM, phylo_nest_1trait = phylo.beta.SNE, phylo_pnest_1trait = pnest),
  func.sor_pam1trait[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_1trait= phylo.beta.SOR, func_turn_1trait = phylo.beta.SIM, func_nest_1trait = phylo.beta.SNE, func_pnest_1trait = pnest))

# at least five traits
beta_5trait <- bind_cols(
  phylo.sor_pam5trait[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_5trait= phylo.beta.SOR, phylo_turn_5trait = phylo.beta.SIM, phylo_nest_5trait = phylo.beta.SNE, phylo_pnest_5trait = pnest),
  func.sor_pam5trait[, c( 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_5trait= phylo.beta.SOR, func_turn_5trait = phylo.beta.SIM, func_nest_5trait = phylo.beta.SNE, func_pnest_5trait = pnest))

# all species
beta_alpha6 <- bind_cols(
  phylo.sor[, c(1, 7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(phylo_total_alpha6= phylo.beta.SOR, phylo_turn_alpha6 = phylo.beta.SIM, phylo_nest_alpha6 = phylo.beta.SNE, phylo_pnest_alpha6 = pnest),
  func.sor[, c(7, 5, 6, 12)] %>% 
    as_tibble() %>%
    rename(func_total_alpha6= phylo.beta.SOR, func_turn_alpha6 = phylo.beta.SIM, func_nest_alpha6 = phylo.beta.SNE, func_pnest_alpha6 = pnest))


## mapping spatial patterns 

# the names of beta variables were wrote in the title
titles <- c("A Phylogenetic total beta-diversity", "B Functional total beta-diversity",
            "C Phylogenetic turnover", "D Functional turnover",
            "E Phylogenetic nestedness", "F Functional nestedness",
            "G Phylogenetic nestedness proportion", "H Functional nestedness proportion")


###### beta diversity using species at least one trait
betaobs_1trait <- beta_1trait[ , c(1,2,6,3,7,4,8,5,9)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(betaobs_1trait$phylo_total_1trait, betaobs_1trait$func_total_1trait)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(betaobs_1trait$phylo_turn_1trait, betaobs_1trait$func_turn_1trait)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(betaobs_1trait$phylo_nest_1trait, betaobs_1trait$func_nest_1trait)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(betaobs_1trait$phylo_pnest_1trait, betaobs_1trait$func_pnest_1trait)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)


# generate maps
beta.maps <- list()
for(i in c(1:8)){
  if(i %in% 1:2)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_1trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 3:4)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_1trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 5:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_1trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:8)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_1trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=2, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_Observed_beta_1trait_maps.png", unit="mm", width=120, height=100)



###### beta diversity using species at least five traits
betaobs_5trait <- beta_5trait[ , c(1,2,6,3,7,4,8,5,9)] %>%
  as.data.frame()

# use the same breaks for three facet of beta diversity, with 0.5 as the median of color palette
# the total beta
total <- c(betaobs_5trait$phylo_total_5trait, betaobs_5trait$func_total_5trait)
breaks.total <- round(seq(min(total), max(total), length=5), 2)
# the turnover
turn <- c(betaobs_5trait$phylo_turn_5trait, betaobs_5trait$func_turn_5trait)
breaks.turn <- round(seq(min(turn), max(turn), length=5), 2)
# the nestedness
nestedness <- c(betaobs_5trait$phylo_nest_5trait, betaobs_5trait$func_nest_5trait)
breaks.nest <- round(seq(min(nestedness), max(nestedness), length=5), 2)
# for proportion of beta contributed by nestedness, use 0.5 as the median of color palette
pnest <- c(betaobs_5trait$phylo_pnest_5trait, betaobs_5trait$func_pnest_5trait)
breaks.pnest <- round(c(seq(min(pnest), 0.5, length =3), seq(0.5, max(pnest), length =3)[-1]), 2)


# generate maps
beta.maps <- list()
for(i in c(1:8)){
  if(i %in% 1:2)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_5trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.total, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 3:4)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_5trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.turn, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 5:6)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_5trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.nest, break.type="linear", mymidpoint="auto", 
                       digits=2, limit=0, cols=NULL, style="cont", inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  if(i %in% 7:8)
    beta.map <- mytmap(myraster=raster_cells, myvalue=betaobs_5trait[,c(1, i+1)], polygons=land_bhm, breaks=breaks.pnest, break.type="linear", mymidpoint=0.5, 
                       digits=2, limit=0, cols=NULL, style="cont",inner.margins=c(0.01,0,0.17,0), 
                       title=titles[i], title.size=0.7, title.position = c(0.10, 0.92),
                       legend.position=c(0.05, 0.01), legend.height=-0.6, legend.text.size=0.6)
  beta.maps[[i]] <- beta.map
}

# save figure
beta.maps <- tmap_arrange(beta.maps, nrow=4, ncol=2, outer.margins=0)
tmap_save(beta.maps, file="results/Fig.S_Observed_beta_5trait_maps.png", unit="mm", width=120, height=100)




###################
# figure showing correlations between beta variables calculated with all species and those species with at least 1, 3 , 5 traits

# set data.frames to store results
components <- c("Nestedness proportion","Nestedness","Turnover","Total")
dimensions <- c("Phylogenetic","Functional")
beta.cor.all.1trait <- expand.grid(dimensions,components)
beta.cor.all.5trait <- expand.grid(dimensions,components)
id <- c(5,9, 4,8, 3,7, 2,6)

# pearson's correlation between all species and species with at least one trait
match.rows <- intersect(beta_alpha6$ID, beta_1trait$ID)
rows1 <- match(match.rows, beta_alpha6$ID)
rows2 <- match(match.rows, beta_1trait$ID)
for(i in 1:8){
  beta.cor.all.1trait[i,"cor"] <- round(cor(beta_alpha6[rows1, id[i]], beta_1trait[rows2, id[i]]), 2)
}

# Pearson's correlation between all species and species with at least five traits
match.rows <- intersect(beta_alpha6$ID, beta_5trait$ID)
rows1 <- match(match.rows, beta_alpha6$ID)
rows2 <- match(match.rows, beta_5trait$ID)
for(i in 1:8){
  beta.cor.all.5trait[i,"cor"] <- round(cor(beta_alpha6[rows1, id[i]], beta_5trait[rows2, id[i]]), 2)
}

# combine different comparison
beta.cor.all.xtrait <- bind_rows(beta.cor.all.1trait %>% mutate(panel = "All species vs >= 1 traits"),
                                 beta.cor.all.5trait %>% mutate(panel = "All species vs >= 5 traits")) %>%
  mutate(panel = factor(panel, levels = c("All species vs >= 1 traits", "All species vs >= 5 traits")))

# generate figure
library(egg) # a package used to generate tags
plot_beta_cor <- ggplot(data=beta.cor.all.xtrait, aes(Var1, Var2, fill= cor)) +
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

ggsave(file="results/Fig.S_correlation_beta_allSpecies.vs.species1or5traits.png", unit="mm",width=100, height=70)



###########################
## Scatter plot comparing functional and phylogenetic beta diversity based on all species and species with at least one and five traits

## organize beta diversity variable using species with at least one or five traits, or all species
# at least one trait
phylo_func_beta_1trait <- beta_1trait %>%
  pivot_longer(col = 2:9, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("_1trait", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  pivot_wider(names_from = beta_facet, values_from = beta_values) %>%
  mutate(traits = "Species with >= 1 trait")

# at least five traits
phylo_func_beta_5trait <- beta_5trait %>%
  pivot_longer(col = 2:9, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("_5trait", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  pivot_wider(names_from = beta_facet, values_from = beta_values) %>%
  mutate(traits = "Species with >= 5 trait")

# all species
phylo_func_beta_all <- beta_alpha6 %>%
  pivot_longer(col = 2:9, names_to = "beta_metric", values_to = "beta_values") %>%
  mutate(beta_metric = gsub("_alpha6", "", beta_metric)) %>%
  separate(beta_metric, into = c("beta_facet", "beta_comp")) %>%
  pivot_wider(names_from = beta_facet, values_from = beta_values) %>%
  mutate(traits = "All species")


phylo_func_beta <- bind_rows(phylo_func_beta_all, phylo_func_beta_1trait, phylo_func_beta_5trait) %>%
  mutate(beta_comp = case_when(beta_comp == "total" ~ "Total beta-diversity",
                               beta_comp == "turn" ~ "Turnover",
                               beta_comp == "nest" ~ "Nestedness",
                               beta_comp == "pnest" ~ "Nestedness proportion"),
         beta_comp = factor(beta_comp, levels = c("Total beta-diversity", "Turnover", "Nestedness", "Nestedness proportion"))) %>%
  left_join(beta_env %>% dplyr::select(ID, longitude, latitude))


# calculate correlation between phylognetic and functional beta-diversity and significance using modified t-test
betacomp_traits <- phylo_func_beta %>% distinct(beta_comp, traits) %>%
  mutate(cor.r = NA,
         cor.p = NA)

for(i in 1:nrow(betacomp_traits)){
  phylo_func_beta_sub <- phylo_func_beta %>% dplyr::filter(beta_comp == betacomp_traits$beta_comp[i] & traits == betacomp_traits$traits[i])
  xy.cor <- modified.ttest(x=phylo_func_beta_sub$phylo, y=phylo_func_beta_sub$func, coords=phylo_func_beta_sub[, c("longitude", "latitude")], nclass=NULL)
  betacomp_traits$cor.r[i] = xy.cor$corr
  betacomp_traits$cor.p[i] = xy.cor$p.value
}

# prepare correlation text added into figures
betacomp_traits <- betacomp_traits %>%
  mutate(r_text = ifelse(cor.p<= 0.001, paste0('italic(r) == ', round(cor.r, 3), '*"***"'), 
                         ifelse(cor.p<= 0.01 & cor.p> 0.001, paste0('italic(r) == ', round(cor.r, 3), '*"**"'),
                                ifelse(cor.p<= 0.05 & cor.p> 0.01, paste0('italic(r) == ', round(cor.r, 3), '*"*"'),
                                       paste0('italic(r) == ', round(r2, 3), "^ns")))))

# generate scatter figures 
ggplot(phylo_func_beta) +
  facet_grid(traits ~ beta_comp, scales = "fixed") +
  geom_point(aes(phylo, func), alpha = 0.2, size = 0.4) +
  geom_smooth(aes(phylo, func), method = "lm", se = FALSE, color = "red", size = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = "blue", lty = 2, size = 0.5) +
  geom_text(data = betacomp_traits, 
            aes(0.4, Inf, label = r_text), hjust = 0.5, vjust = 2.5, size = 2.81, parse =TRUE) + 
  coord_fixed() +
  scale_x_continuous(limits = c(0, 0.95), breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_y_continuous(limits = c(0, 0.95), breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(x = "Phylogenetic", y = "Functional") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8, colour = "black"),
        strip.text = element_text(size = 9))

ggsave(file="results/Fig.S_compare_phylo_func_beta_traits.png.png", unit="mm",width=180, height=140, dpi = 600)



############################
# phylogenetic signals

# remove some gymnosperms in the trait data
trait_obs <- trait_obs %>% dplyr::filter(species %in% tree$tip.label)

################
# calculate phylogenetic signals of imputed traits using the same number of species with observed traits
# sample 200 time and calculate the mean values

# A function to sample imputed functional traits and calculate phylogenetic signals
phylosig_sample <- function(trait_obs = trait_obs, trait_pre = trait, tree = tree){
  trait_name <- unique(trait_obs$trait_name)
  phylosig_pre <- data.frame(trait_name = trait_name, 
                             lambda = NA, lambda_P = NA, K = NA, K_P = NA)
  for(i in 1:length(trait_name)){
    # the imputed trait values
    trait_pre_sub <- trait_pre %>% pull(trait_name[i])
    trait_values <- setNames(trait_pre_sub, trait_pre$species)
    # sample same number of species with observed traits from the imputed traits
    nspp_obs <- trait_obs %>% dplyr::filter(trait_name == trait_name[i]) %>% nrow()
    trait_values_sample <- sample(trait_values, nspp_obs)
    # calculate phylogenetic signals
    lambda_trait <- phylosig(tree, trait_values_sample, method="lambda", test=TRUE)
    K_trait <- phylosig(tree, trait_values_sample, method="K", test=TRUE)
    phylosig_pre[i, "lambda"] <- lambda_trait$lambda
    phylosig_pre[i, "lambda_P"] <- lambda_trait$P
    phylosig_pre[i, "K"] <- K_trait$K
    phylosig_pre[i, "K_P"] <- K_trait$P
    print(i)
  }
  return(phylosig_pre)
}

# set parallel computing
print("n_cores")
print(detectCores()) #the number of cores
no_cores <- 25
cl <- makeCluster(no_cores, type="PSOCK") # Initiate cluster
clusterExport(cl, varlist = c("trait_obs", "trait", "tree", "phylosig_sample"))
clusterEvalQ(cl, c(library(tidyverse), library(phytools)))

# calculate phylogenetic signals
get_phylosig_tempo <- function(i) phylosig_sample(trait_obs = trait_obs, trait_pre = trait, tree = tree)
phylosig_pre_trait <- parLapply(cl, 1:200, get_phylosig_tempo)

save(phylosig_pre_trait, file = "results/trait_phylosig_obs_pre.RDATA")
stopCluster(cl)


# calculate the phylogenetic signal of traits before imputation
# a data frame to save phylogenetic signals
trait_phylosig <- data.frame(trait_name = unique(trait_obs$trait_name), 
                             lambda = NA, lambda_P = NA, K = NA, K_P = NA)

for(i in 1:8){
  trait_obs_sub <- trait_obs %>% dplyr::filter(trait_name == trait_phylosig[i, 1])
  trait_obs_values <- setNames(trait_obs_sub$trait_values, trait_obs_sub$species)
  lambda_trait <- phylosig(tree, trait_obs_values, method="lambda", test=TRUE)
  K_trait <- phylosig(tree, trait_obs_values, method="K", test=TRUE)
  trait_phylosig[i, "lambda"] <- lambda_trait$lambda
  trait_phylosig[i, "lambda_P"] <- lambda_trait$P
  trait_phylosig[i, "K"] <- K_trait$K
  trait_phylosig[i, "K_P"] <- K_trait$P
  print(i)
}

save(trait_phylosig, phylosig_pre_trait, file = "results/trait_phylosig_obs_pre.RDATA")


## calculate mean and sd of phylogenetic signals using traits after imputation
phylosig_pre_trait <- bind_rows(phylosig_pre_trait) %>%
  as_tibble()

phylosig_pre_trait_summary <- phylosig_pre_trait %>%
  group_by(trait_name) %>%
  summarise(lambda_mean = mean(lambda),
            lambda_sd = sd(lambda),
            n_lambda_sig = sum(lambda_P<0.001),
            K_mean = mean(K),
            K_sd = sd(K),
            n_K_sig = sum(K_P<0.05))

# format phylogenetic signals.
# Choose Pagel's lambda as this metric seems more robust to polytomis (Molina-Venegas and Rodriguez, 2017).
phylosig_obs_pre_trait <- phylosig_pre_trait_summary %>%
  dplyr::select(trait_name, lambda_mean, lambda_sd, n_lambda_sig) %>%
  mutate(lambda_pre_mean_sd = paste0(round(lambda_mean, 3), " (", round(lambda_sd, 3) ,")")) %>%
  left_join(trait_phylosig %>% dplyr::select(trait_name, lambda, lambda_P)) %>%
  dplyr::select(trait_name, lambda_obs = lambda, lambda_P, lambda_pre_mean_sd, n_lambda_sig) %>%
  mutate(lambda_obs = round(lambda_obs, 3))

write.csv(phylosig_obs_pre_trait, file = "results/Table.S_trait_phylosig.csv", row.names = FALSE)
