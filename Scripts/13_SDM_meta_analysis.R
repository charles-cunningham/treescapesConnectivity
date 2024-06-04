# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: SDM meta analysis
#
# Script Description: Analysis of SDM inlabru model outputs using brms package

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Run once for R 4.2.2 to get the brms package working
# Install RTools 42
# remove.packages(c("StanHeaders", "rstan", "brms"))
# if (file.exists(".RData")) file.remove(".RData")
# RESTART R
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# example(stan_model, package = "rstan", run.dontrun = TRUE) # This checks rstan and the C++ compiler are correctly installed
# install.packages("brms")
# RESTART R

# Load packages
library(brms)
library(terra)
library(tidyverse)
library(tidybayes)
library(ggridges)

# SET PARAMETERS ------------------------------------

# DATA FILES ------------------------------------------

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# DESCRIPTIVE TAXA PLOTS ----------------------------------

# Box plot of connectivity model effect sizes

# Plot
connEffPlot <- meta_df %>%
  filter(abs(mean_connectivity) < 2.5) %>% # lots of outliers
  ggplot() +
  geom_boxplot(aes(x = taxa, y = mean_connectivity)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip()

# Save
ggsave(filename = paste0("../Writing/Plots/", "Raw_connectivity_effect_sizes.png"),
       connEffPlot,
       dpi = 600,
       units = "px", width = 6000, height = 5000)

### META ANALYSIS --------------------------------------
# In this stage we want to check overall data patterns,
# partitioning of variance between taxa and species levels etc.
# Run all brm functions first in a single section. It takes a while
# so can save the outputs so only need to run once.

# RUN BAYESIAN MODELS ---------------------------------

# Fit Bayesian meta analysis model
# N.B. Species effect nested within taxa effect.
# (i.e. taxa + taxa:species random effects)

### Connectivity

# All species
connAll_brms <- brm(data = meta_df,
                    family = gaussian,
                    mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa ),
                    prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 0.5), class = sd)),
                    iter = 10000,
                    warmup = 5000, 
                    control=list(adapt_delta = 0.99,
                                 stepsize = 0.1,
                                 max_treedepth = 15),
                    cores = 4,
                    chains = 4)

# Broadleaf species only
connBF_brms <- brm(data = meta_df[meta_df$broadleafAssociation == "Y", ],
                   family = gaussian,
                   mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa),
                    prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 0.5), class = sd)),
                   iter = 10000,
                   warmup = 5000, 
                   control=list(adapt_delta = 0.99,
                                stepsize = 0.1,
                                max_treedepth = 15),
                   cores = 4,
                   chains = 4)

# Coniferous species only
connCF_brms <- brm(data = meta_df[meta_df$coniferousAssociation == "Y", ],
                   family = gaussian,
                   mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa ),
                   prior = c(prior(normal(0, 1), class = Intercept),
                             prior(cauchy(0, 0.5), class = sd)),
                   iter = 10000,
                   warmup = 5000, 
                   control=list(adapt_delta = 0.99,
                                stepsize = 0.1,
                                max_treedepth = 15),
                   cores = 4,
                   chains = 4)

# Woodland avoiding species only
connOpen_brms <- brm(data = meta_df[meta_df$openAssociation == "Y", ],
                     family = gaussian,
                     mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa ),
                     prior = c(prior(normal(0, 1), class = Intercept),
                               prior(cauchy(0, 0.5), class = sd)),
                     iter = 10000,
                     warmup = 5000, 
                     control=list(adapt_delta = 0.99,
                                  stepsize = 0.1,
                                  max_treedepth = 15),
                     cores = 4,
                     chains = 4)

### Cover

# Broadleaf species only
coverBF_brms <- brm(data = meta_df[meta_df$broadleafAssociation == "Y", ],
                    family = gaussian,
                    mean_coverBF | se(sd_coverBF) ~ 1 + (1 | taxa ),
                    prior = c(prior(normal(0, 1), class = Intercept),
                            prior(cauchy(0, 0.5), class = sd)),
                    iter = 10000,
                    warmup = 5000, 
                    control=list(adapt_delta = 0.99,
                                 stepsize = 0.1,
                                 max_treedepth = 15),
                    cores = 4,
                    chains = 4)

# Coniferous species only
coverCF_brms <- brm(data = meta_df[meta_df$coniferousAssociation == "Y", ],
                     family = gaussian,
                     mean_coverCF | se(sd_coverCF) ~ 1 + (1 | taxa ),
                     prior = c(prior(normal(0, 1), class = Intercept),
                               prior(cauchy(0, 0.5), class = sd)),
                    iter = 10000,
                    warmup = 5000, 
                    control=list(adapt_delta = 0.99,
                                 stepsize = 0.1,
                                 max_treedepth = 15),
                    cores = 4,
                    chains = 4)

### Cover:connectivity interaction

# Woodland species only
intBF_brms <- brm(data = meta_df[meta_df$broadleafAssociation == "Y", ],
                  family = gaussian,
                  mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(cauchy(0, 0.5), class = sd)),
                  iter = 10000,
                  warmup = 5000, 
                  control=list(adapt_delta = 0.99,
                               stepsize = 0.1,
                               max_treedepth = 15),
                  cores = 4,
                  chains = 4)
# Woodland species only
intCF_brms <- brm(data = meta_df[meta_df$coniferousAssociation == "Y", ],
                        family = gaussian,
                        mean_CFconnINT | se(sd_CFconnINT) ~ 1 + (1 | taxa ),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 0.5), class = sd)),
                  iter = 10000,
                  warmup = 5000,
                  control=list(adapt_delta = 0.99,
                               stepsize = 0.1,
                               max_treedepth = 15),
                  cores = 4,
                  chains = 4)

### Save

# List of brms objects
brmsList <- c("connAll_brms", "connBF_brms", "connCF_brms", "connOpen_brms",
              "coverBF_brms", "coverCF_brms",
              "intBF_brms", "intCF_brms")

# Save brms objects
save(list = brmsList,
     file = "../Data/Species_data/brms_objects.RData")