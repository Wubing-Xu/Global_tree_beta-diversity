##################################################################################################################
## Organize observed and expected three facets of beta-diversity within regions of 3*3 grid-cells (8 nearest neighboring cells around the focal),
# and then calculate the deviation between the observed and expected phylogenetic and functional beta-diversity
##################################################################################################################

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/AU/global_tree_beta_2022",
                  "IDIVTS01" = "H:/wubing/AU/global_tree_beta_2022")
setwd(path2wd)


# load packages
needed_libs <- c("tidyverse", "betapart", "picante", "ape", "phytools", "raster", "spdep", "parallel")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("intermediate_results/observed_beta_sorensen.RDATA")
load("intermediate_results/random_beta_phylo_func_nn8.RDATA")



# combine observed species, phylogenetic and functional beta diversity
spesor_nn8 <- spe.sor_nn8 %>% 
  as_tibble() %>%
  mutate(alpha.max = shared + max.notshared,
         alpha.min = shared + min.notshared,
         diff.alpha = max.notshared - min.notshared) %>%
  rename(beta.turn = beta.SIM, beta.nest = beta.SNE, beta.total= beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("spe", names(.)[-c(1:4)], sep = ".")))

phylosor_nn8 <- phylo.sor_nn8 %>% 
  as_tibble() %>%
  mutate(alpha.max = shared + max.notshared,
         alpha.min = shared + min.notshared,
         diff.alpha = max.notshared - min.notshared) %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("phylo", names(.)[-c(1:4)], sep = ".")))

funcsor_nn8 <- func.sor_nn8 %>% 
  as_tibble() %>%
  mutate(alpha.max = shared + max.notshared,
         alpha.min = shared + min.notshared,
         diff.alpha = max.notshared - min.notshared) %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("func", names(.)[-c(1:4)], sep = ".")))

betasor_nn8 <- bind_cols(spesor_nn8, phylosor_nn8[-c(1:4)], funcsor_nn8[-c(1:4)])


## organize mean and sd of expected phylogenetic and functional beta from null model
phylosor_nn8_null_mean <- phylo.sor_nn8_null_mean %>% 
  as_tibble() %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("phylo", names(.)[-c(1:4)], sep = ".")))

phylosor_nn8_null_sd <- phylo.sor_nn8_null_sd %>% 
  as_tibble() %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("phylo", names(.)[-c(1:4)], sep = ".")))

funcsor_nn8_null_mean <- func.sor_nn8_null_mean %>% 
  as_tibble() %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("func", names(.)[-c(1:4)], sep = ".")))

funcsor_nn8_null_sd <- func.sor_nn8_null_sd %>% 
  as_tibble() %>%
  rename(beta.turn = phylo.beta.SIM, beta.nest = phylo.beta.SNE, beta.total= phylo.beta.SOR) %>%
  relocate(c(beta.total, beta.turn, beta.nest, pnest), .after = gamma) %>%
  setNames(c(names(.)[1:4], paste("func", names(.)[-c(1:4)], sep = ".")))


## calculate deviation between the observed and expected beta 
phylosor_nn8_devia <- phylosor_nn8_null_mean
phylosor_nn8_devia[, 5:8] <- phylosor_nn8[, 5:8] - phylosor_nn8_null_mean[, 5:8]
phylosor_nn8_devia[, 9:15] <- (phylosor_nn8[, 9:15] - phylosor_nn8_null_mean[, 9:15])/phylosor_nn8[, 9:15]

funcsor_nn8_devia <- funcsor_nn8_null_mean
funcsor_nn8_devia[, 5:8] <- funcsor_nn8[, 5:8] - funcsor_nn8_null_mean[, 5:8]
funcsor_nn8_devia[, 9:15] <- (funcsor_nn8[, 9:15] - funcsor_nn8_null_mean[, 9:15])/funcsor_nn8[, 9:15]


## calculate SES of beta as the difference between the observed and mean expected beta dividing by the sd of expected beta
phylosor_nn8_ses <- phylosor_nn8_null_mean
phylosor_nn8_ses[, 5:15] <- (phylosor_nn8[, 5:15] - phylosor_nn8_null_mean[, 5:15])/phylosor_nn8_null_sd[, 5:15]

funcsor_nn8_ses <- funcsor_nn8_null_mean
funcsor_nn8_ses[, 5:15] <- (funcsor_nn8[, 5:15] - funcsor_nn8_null_mean[, 5:15])/funcsor_nn8_null_sd[, 5:15]


##Save results
save(spesor_nn8, phylosor_nn8, funcsor_nn8, betasor_nn8, 
     phylosor_nn8_null_mean, phylosor_nn8_null_sd,
     funcsor_nn8_null_mean, funcsor_nn8_null_sd,
     phylosor_nn8_devia, funcsor_nn8_devia,
     phylosor_nn8_ses, funcsor_nn8_ses,
     file="intermediate_results/beta_observed_expected_deviation_nn8.RDATA")
