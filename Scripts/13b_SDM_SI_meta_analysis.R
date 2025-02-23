# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: SDM meta analysis
#
# Script Description: Analysis of SDM inlabru model outputs using brms package

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Run once for R 4.4.2 to get the brms package working (belt and braces)
# Install RTools
# remove.packages(c("StanHeaders", "rstan", "brms"))
# if (file.exists(".RData")) file.remove(".RData")
# RESTART R
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# options(mc.cores = parallel::detectCores())
# example(stan_model, package = "rstan", run.dontrun = TRUE) # This checks rstan and the C++ compiler are correctly installed
# install.packages("brms")
# rstan_options(auto_write = TRUE)
# RESTART R

# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)

# SET PARAMETERS ------------------------------------

# DATA FILES ------------------------------------------

# SDM data folder
dataDir <- paste0("../Data/Species_data/SDMs/")

### META ANALYSIS --------------------------------------

# UNDER 30% COVER --------------------------------------

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_below_30",
                   ".RData"))

# Fit Bayesian meta analysis model
# N.B. Species effect nested within taxa effect
# (i.e. taxa + taxa:species random effects)

### Connectivity

# Broadleaf species only
connBF_brms_under30 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                           family = gaussian,
                           mean_connectivity | se(sd_connectivity) ~
                             1 + (1 | taxa) + (1 | taxa:species),
                           prior = c(prior(normal(0, 1), class = Intercept),
                                     prior(cauchy(0, 1), class = sd)),
                           iter = 10000,
                           warmup = 5000, 
                           control=list(adapt_delta = 0.999,
                                        stepsize = 0.001,
                                        max_treedepth = 20),
                           cores = 4,
                           chains = 4)

### Cover

# Broadleaf species only
coverBF_brms_under30 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                            family = gaussian,
                            mean_coverBF | se(sd_coverBF) ~
                              1 + (1 | taxa) + (1 | taxa:species),
                            prior = c(prior(normal(0, 1), class = Intercept),
                                      prior(cauchy(0, 1), class = sd)),
                            iter = 10000,
                            warmup = 5000, 
                            control=list(adapt_delta = 0.999,
                                         stepsize = 0.001,
                                         max_treedepth = 20),
                            cores = 4,
                            chains = 4)

### Cover:connectivity interaction

# Woodland species only
intBF_brms_under30 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                          family = gaussian,
                          mean_BFconnINT | se(sd_BFconnINT) ~
                            1 + (1 | taxa) + (1 | taxa:species),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          iter = 10000,
                          warmup = 5000, 
                          control=list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20),
                          cores = 4,
                          chains = 4)

### Save

# List of brms objects
brmsList <- c("connBF_brms_under30", "coverBF_brms_under30", "intBF_brms_under30")

# Save brms objects
save(list = brmsList,
     file = "../Data/Species_data/brms_objects_below_30.RData")

# QUADRATIC ANALYSIS --------------------------------------------

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_quadratic",
                   ".RData"))

# Fit Bayesian meta analysis model
# N.B. Species effect nested within taxa effect
# (i.e. taxa + taxa:species random effects)

### Connectivity

# Broadleaf species only
connBF_brms_quad <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                        family = gaussian,
                        mean_connectivity | se(sd_connectivity) ~
                          1 + (1 | taxa) + (1 | taxa:species),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 10000,
                        warmup = 5000, 
                        control=list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                        cores = 4,
                        chains = 4)

### Cover

# Broadleaf species only
coverBF_brms_quad <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                         family = gaussian,
                         mean_coverBF | se(sd_coverBF) ~
                           1 + (1 | taxa) + (1 | taxa:species),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                         iter = 10000,
                         warmup = 5000, 
                         control=list(adapt_delta = 0.999,
                                      stepsize = 0.001,
                                      max_treedepth = 20),
                         cores = 4,
                         chains = 4)

# Woodland species only
quadBF_brms_quad <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                        family = gaussian,
                        mean_coverBF_quad | se(sd_coverBF_quad ) ~
                          1 + (1 | taxa) + (1 | taxa:species),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 10000,
                        warmup = 5000, 
                        control=list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                        cores = 4,
                        chains = 4)

### Cover:connectivity interaction

# Woodland species only
intBF_brms_quad <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                       family = gaussian,
                       mean_BFconnINT | se(sd_BFconnINT) ~
                         1 + (1 | taxa) + (1 | taxa:species),
                       prior = c(prior(normal(0, 1), class = Intercept),
                                 prior(cauchy(0, 1), class = sd)),
                       iter = 10000,
                       warmup = 5000, 
                       control=list(adapt_delta = 0.999,
                                    stepsize = 0.001,
                                    max_treedepth = 20),
                       cores = 4,
                       chains = 4)

### Save

# List of brms objects
brmsList <- c("connBF_brms_quad", "coverBF_brms_quad",
              "intBF_brms_quad", "quadBF_brms_quad")

# Save brms objects
save(list = brmsList,
     file = "../Data/Species_data/brms_objects_quadratic.RData")

# QUARTILE ANALYSIS --------------------------------------------

# 75-100% QUARTILE

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_quartile_75-100",
                   ".RData"))

# Save species to use for other quadratic runs
# N.B. Standardise to these as if species appears here, it appears for all
quadSpecies <- unique(SI_df$species)

# Connectivity
connBF_brms_75_100 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                          family = gaussian,
                          mean_connectivity | se(sd_connectivity) ~
                            1 + (1 | taxa) + (1 | taxa:species),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          iter = 10000,
                          warmup = 5000, 
                          control=list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20),
                          cores = 4,
                          chains = 4)

# 0-25% QUARTILE

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_quartile_0-25",
                   ".RData"))

# Filter species to quadSpecies
SI_df <- SI_df %>%
  filter(species %in% quadSpecies)

# Connectivity
connBF_brms_0_25 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                        family = gaussian,
                        mean_connectivity | se(sd_connectivity) ~
                          1 + (1 | taxa) + (1 | taxa:species),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        iter = 10000,
                        warmup = 5000, 
                        control=list(adapt_delta = 0.999,
                                     stepsize = 0.001,
                                     max_treedepth = 20),
                        cores = 4,
                        chains = 4)

# 25-50% QUARTILE

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_quartile_25-50",
                   ".RData"))

# Filter species
SI_df <- SI_df %>%
  filter(species %in% quadSpecies)

# Connectivity
connBF_brms_25_50 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                   family = gaussian,
                   mean_connectivity | se(sd_connectivity) ~
                     1 + (1 | taxa) + (1 | taxa:species),
                   prior = c(prior(normal(0, 1), class = Intercept),
                             prior(cauchy(0, 1), class = sd)),
                   iter = 10000,
                   warmup = 5000, 
                   control=list(adapt_delta = 0.999,
                                stepsize = 0.001,
                                max_treedepth = 20),
                   cores = 4,
                   chains = 4)

# 50-75% QUARTILE

# Load data
load(file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   "SI_quartile_50-75",
                   ".RData"))

# Filter species
SI_df <- SI_df %>%
  filter(species %in% quadSpecies)

# Connectivity
connBF_brms_50_75 <- brm(data = SI_df[SI_df$broadleafAssociation == "Y", ],
                   family = gaussian,
                   mean_connectivity | se(sd_connectivity) ~
                     1 + (1 | taxa) + (1 | taxa:species),
                   prior = c(prior(normal(0, 1), class = Intercept),
                             prior(cauchy(0, 1), class = sd)),
                   iter = 10000,
                   warmup = 5000, 
                   control=list(adapt_delta = 0.999,
                                stepsize = 0.001,
                                max_treedepth = 20),
                   cores = 4,
                   chains = 4)

# List of brms objects
brmsList <- c("connBF_brms_75_100", "connBF_brms_0_25",
              "connBF_brms_25_50", "connBF_brms_50_75")

# Save brms objects
save(list = brmsList,
     file = "../Data/Species_data/brms_objects_quart.RData")
