# HEADER --------------------------------------------
#
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

# Categorisation of effect sizes for woodland and all species

# Summarise effect sizes of all species SDMs by BF woodland species
groupedSummary <- meta_df %>%
  # Group by taxa group and effect category...
  group_by(taxa, connectivity_sig, woodSp ) %>% # Change to CFwoodSp to see coniferous associations
  # ... and count total number
  summarise(nuSpecies = length(species)) %>%
  # Ungroup for ...
  ungroup %>%
  # ... proportion
  mutate(freq = prop.table(nuSpecies), .by = c(taxa, woodSp))

# Plot
ggplot(groupedSummary) +
  geom_bar(aes(x = taxa, y = nuSpecies, fill = connectivity_sig),
           position = "stack",
           stat = "identity") +
  facet_wrap( ~ woodSp) + 
coord_flip() 

ggplot(groupedSummary) +
  geom_bar(aes(x = taxa, y = freq, fill = connectivity_sig),
           position = "stack",
           stat = "identity") +
  facet_wrap( ~ woodSp) +
  coord_flip() 

### META ANALYSIS --------------------------------------
# In this stage we want to check overall data patterns,
# partitioning of variance between taxa and species levels etc.
# Run all brm functions first in a single section. It takes a while
# so can save the outputs so only need to run once.

# RUN BAYESIAN MODELS ---------------------------------

# # Fit Bayesian meta analysis model
# # N.B. Species effect nested within taxa effect.
# # (i.e. taxa + taxa:species random effects)
# 
# ### Connectivity
# 
# # All species
# connAll_brms <- brm(data = meta_df,
#               family = gaussian,
#               mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
#               prior = c(prior(normal(0, 1), class = Intercept),
#                         prior(cauchy(0, 1), class = sd)),
#               iter = 4000,
#               cores = 4,
#               chains = 4)
# # Woodland species only
# connWood_brms <- brm(data = meta_df[meta_df$woodSp == "Prefer woodland", ],
#                  family = gaussian,
#                  mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
#                  prior = c(prior(normal(0, 1), class = Intercept),
#                            prior(cauchy(0, 1), class = sd)),
#                  iter = 4000,
#                  cores = 4,
#                  chains = 4)
# # Woodland avoiding species only
# connNotWood_brms <- brm(data = meta_df[meta_df$woodSp == "Avoid woodland", ],
#                   family = gaussian,
#                   mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
#                   prior = c(prior(normal(0, 1), class = Intercept),
#                             prior(cauchy(0, 1), class = sd)),
#                   iter = 4000,
#                   cores = 4,
#                   chains = 4)
# 
# ### Cover
# 
# # All species
# coverAll_brms <- brm(data = meta_df,
#                  family = gaussian,
#                  mean_coverBF | se(sd_coverBF) ~ 1 + (1 | taxa / species),
#                  prior = c(prior(normal(0, 1), class = Intercept),
#                            prior(cauchy(0, 1), class = sd)),
#                  iter = 4000,
#                  cores = 4,
#                  chains = 4)
# # Woodland species only
# coverWood_brms <- brm(data = meta_df[meta_df$woodSp == "Prefer woodland", ],
#                   family = gaussian,
#                   mean_coverBF | se(sd_coverBF) ~ 1 + (1 | taxa / species),
#                   prior = c(prior(normal(0, 1), class = Intercept),
#                             prior(cauchy(0, 1), class = sd)),
#                   iter = 4000,
#                   cores = 4,
#                   chains = 4)
# # Woodland species only
# coverNotWood_brms <- brm(data = meta_df[meta_df$woodSp == "Avoid woodland", ],
#                      family = gaussian,
#                      mean_coverBF | se(sd_coverBF) ~ 1 + (1 | taxa / species),
#                      prior = c(prior(normal(0, 1), class = Intercept),
#                                prior(cauchy(0, 1), class = sd)),
#                      iter = 4000,
#                      cores = 4,
#                      chains = 4)
# 
# ### Cover:connectivity interaction
# 
# # All species
# intAll_brms <- brm(data = meta_df,
#                     family = gaussian,
#                     mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
#                     prior = c(prior(normal(0, 1), class = Intercept),
#                               prior(cauchy(0, 1), class = sd)),
#                     iter = 4000,
#                     cores = 4,
#                     chains = 4)
# # Woodland species only
# intWood_brms <- brm(data = meta_df[meta_df$woodSp == "Prefer woodland", ],
#                      family = gaussian,
#                      mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
#                      prior = c(prior(normal(0, 1), class = Intercept),
#                                prior(cauchy(0, 1), class = sd)),
#                      iter = 4000,
#                      cores = 4,
#                      chains = 4)
# # Woodland species only
# intNotWood_brms <- brm(data = meta_df[meta_df$woodSp == "Avoid woodland", ],
#                         family = gaussian,
#                         mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
#                         prior = c(prior(normal(0, 1), class = Intercept),
#                                   prior(cauchy(0, 1), class = sd)),
#                         iter = 4000,
#                         cores = 4,
#                         chains = 4)
# 
# ### Save
# 
# # List of brms objects
# brmsList <- c("connAll_brms", "connWood_brms", "connNotWood_brms",
#               "coverAll_brms", "coverWood_brms", "coverNotWood_brms",
#               "intAll_brms", "intWood_brms", "intNotWood_brms")
# 
# # Save brms objects
# save(list = brmsList,
#      file = "../Data/Species_data/brms_objects.RData")

# ANALYSIS ---------------------------------------------------------------

# Load brms objects
load(file = "../Data/Species_data/brms_objects.RData")

# Print summary for all species and woodland species
summary(coverAll_brms)


# Pairs plot
#pairs(metaAn)

# Posterior check
#pp_check(metaAnWood)

# Random effects
#ranef(metaAn)

### Plots

# Plot tau for taxa and species
# (variance between species larger and more consistent than taxa)

# Set plot titles
titles <- c( "All", "Wood", "Not wood")

# Loop through deffirent meta analyses
for (i in 1:length(metaAnList)) {
  
  # Extract draws, subset to sd, and process for plot
  tauPlot <- as_draws_df(metaAnList[[i]]) %>%
  select(starts_with("sd")) %>%
  gather(key, tau) %>%
  mutate(key = str_remove(key, "sd_") %>% str_remove(., "__Intercept")) %>%
  
  # Plot
  ggplot(aes(x = tau, fill = key)) +
  geom_density(color = "transparent", alpha = 2 / 3) +
  scale_fill_viridis_d(NULL, end = .85) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(expression(tau)) +
  ggtitle(titles[i]) +
  theme(panel.grid = element_blank())

# Save
ggsave(filename = paste0("../Writing/Plots/", "Meta_analysis", titles[i], "_tau.png"),
       tauPlot,
       dpi = 600,
       units = "px", width = 6000, height = 5000)

}

# Plot intercepts

# Loop through meta analysis subsets
for (i in 1:length(metaAnList)) {

  # Assign meta analysis i
  metaAn <- connWood_brms #metaAnList[[i]]

  # Extract draws grouping by taxa
  study.draws <- spread_draws(metaAn, r_taxa[taxa, ], b_Intercept) %>%
    mutate(b_Intercept = r_taxa + b_Intercept)
  
  # Extract draws (all pooled)
  pooled.effect.draws <- spread_draws(metaAn, b_Intercept) %>%
    mutate(taxa = "Pooled Effect")
  
  # Bind taxa draws and pooled draws together
  forest.data <- bind_rows(study.draws,
                           pooled.effect.draws) %>%
    ungroup() %>%
    mutate(taxa = str_replace_all(taxa, "[.]", " ")) %>%
    mutate(taxa = reorder(taxa, b_Intercept))
  
  # Group by taxa and summarize by intercept
  forest.data.summary <- group_by(forest.data, taxa) %>%
    mean_qi(b_Intercept)
  
  # Plot
  taxaSummaries <- ggplot(data = forest.data,
                          aes(b_Intercept,
                              relevel(taxa, "Pooled Effect",
                                      after = Inf),
                              fill = taxa)) +
    
    # Add vertical lines for pooled effect and CI
    geom_vline(xintercept = fixef(metaAn)[1, 1],
               color = "grey",
               size = 1) +
    geom_vline( xintercept = fixef(metaAn)[1, 3:4],
                color = "grey",
                linetype = 2) +
    geom_vline(xintercept = 0,
               color = "black",
               size = 1) +
    
    # Add densities
    geom_density_ridges(scale = 0.95,
                        rel_min_height = 0.001) +
    geom_pointinterval(data = forest.data.summary,
                       size = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Add labels
    labs(x = "Mean effect of connectivity", # summary measure
         y = element_blank()) +
    xlim(-0.7, 0.9) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20))
  
  # Save
  ggsave(
    filename = paste0("../Writing/Plots/", "Meta_analysis",
                      titles[i], "_taxaSummaries.png"),
    taxaSummaries,
    dpi = 600,
    units = "px",
    width = 4000,
    height = 5000
  )
}

#??? freshwater species bottom 3!
hypothesis(metaAnList[[3]], "Intercept > 0")



#####




# Print summary for all species and woodland species
summary(metaAnIntAll)
summary(metaAnIntWood)
summary(metaAnIntNotWood)

# Pairs plot
#pairs(metaAn)

# Posterior check
#pp_check(metaAnWood)

# Random effects
#ranef(metaAn)

### Plots

# Plot tau for taxa and species
# (variance between species larger and more consistent than taxa)

# Set plot titles
titles <- c( "All", "Wood", "Not wood")

# Loop through deffirent meta analyses
for (i in 1:length(metaAnIntList)) {
  
  # Extract draws, subset to sd, and process for plot
  tauPlot <- as_draws_df(metaAnIntList[[i]]) %>%
    select(starts_with("sd")) %>%
    gather(key, tau) %>%
    mutate(key = str_remove(key, "sd_") %>% str_remove(., "__Intercept")) %>%
    
    # Plot
    ggplot(aes(x = tau, fill = key)) +
    geom_density(color = "transparent", alpha = 2 / 3) +
    scale_fill_viridis_d(NULL, end = .85) +
    scale_y_continuous(NULL, breaks = NULL) +
    xlab(expression(tau)) +
    ggtitle(titles[i]) +
    theme(panel.grid = element_blank())
  
  # Save
  ggsave(filename = paste0("../Writing/Plots/", "Meta_analysis", titles[i], "INT_tau.png"),
         tauPlot,
         dpi = 600,
         units = "px", width = 6000, height = 5000)
  
}

# Plot intercepts

# Loop through meta analysis subsets
for (i in 1:length(metaAnIntList)) {
  
  # Assign meta analysis i
  metaAn <- metaAnIntList[[i]]
  
  # Extract draws grouping by taxa
  study.draws <- spread_draws(metaAn, r_taxa[taxa, ], b_Intercept) %>%
    mutate(b_Intercept = r_taxa + b_Intercept)
  
  # Extract draws (all pooled)
  pooled.effect.draws <- spread_draws(metaAn, b_Intercept) %>%
    mutate(taxa = "Pooled Effect")
  
  # Bind taxa draws and pooled draws together
  forest.data <- bind_rows(study.draws,
                           pooled.effect.draws) %>%
    ungroup() %>%
    mutate(taxa = str_replace_all(taxa, "[.]", " ")) %>%
    mutate(taxa = reorder(taxa, b_Intercept))
  
  # Group by taxa and summarize by intercept
  forest.data.summary <- group_by(forest.data, taxa) %>%
    mean_qi(b_Intercept)
  
  # Plot
  taxaSummaries <- ggplot(data = forest.data,
                          aes(b_Intercept,
                              relevel(taxa, "Pooled Effect",
                                      after = Inf),
                              fill = taxa)) +
    
    # Add vertical lines for pooled effect and CI
    geom_vline(xintercept = fixef(metaAn)[1, 1],
               color = "grey",
               size = 1) +
    geom_vline( xintercept = fixef(metaAn)[1, 3:4],
                color = "grey",
                linetype = 2) +
    geom_vline(xintercept = 0,
               color = "black",
               size = 1) +
    
    # Add densities
    geom_density_ridges(scale = 0.95,
                        rel_min_height = 0.001) +
    geom_pointinterval(data = forest.data.summary,
                       size = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Add labels
    labs(x = "Mean effect of connectivity", # summary measure
         y = element_blank()) +
    xlim(-0.7, 0.9) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20))
  
  # Save
  ggsave(
    filename = paste0("../Writing/Plots/", "Meta_analysis",
                      titles[i], "_int_taxaSummaries.png"),
    taxaSummaries,
    dpi = 600,
    units = "px",
    width = 4000,
    height = 5000
  )
}

#??? freshwater species bottom 3!

hypothesis(metaAnList[[3]], "Intercept > 0")


# INTERACTION PLOTS ----------------------------------- IN PROGRESS, temp version for now

# Load scaling parameters
load("../Data/Spatial_data/DataForInlabru/scalingParams.RData")

# How many prediction steps?
nSamp <- 10

### Create unscaled prediction grid

# Load cover/connectivity spatRasters
for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                     i)))
}

# Create dataframes from spatrasters
connW_df <- connW %>%
  as.data.frame(.) %>%
  gather("Year", "Connectivity")
coverW_df <- (coverBF + coverCF) %>%
  as.data.frame(.) %>%
  gather("Year", "Cover")

# Extract max scaled value from connectivity spatraster
maxConnectivity <- connW %>%
  global(., fun = "max", na.rm = TRUE) %>%
  max

# Create unscaled data frame of cover and connectivity values to predict over
# (only broadleaf for now)
BF_pred_df <-  expand.grid(BF_pred = seq(0, 1, by = 1/nSamp),
                           conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))

### Remove prediction grid cells with no actual data points, we want to any predictions to 
### be possible, i.e. high connectivity with very low cover is impossible

# Create data frame of cover and connectivity data points
data_pts <- cbind(coverW_df, connW_df[, "Connectivity"]) %>%
  .[,c(2,3)]
names(pred_pts) <- c("Cover", "Connectivity")

# Convert cover data points to interval bins
coverBins <- findInterval(data_pts[, 1],
                          unique(BF_pred_df[, 1])) %>%
  unique(BF_pred_df[, 1])[.]

# Convert connectivity data points to interval bins
connBins <- findInterval(data_pts[, 2],
                          unique(BF_pred_df[, 2])) %>%
  unique(BF_pred_df[, 2])[.]

# Create a data frame of cross-occurrence of cover and connectivity bins
predBins <- table("Cover" = coverBins,
              "Connectivity" = connBins) %>%
  as.data.frame %>%
  mutate_all(function(x) {as.numeric(as.character(x))} )
  
# Remove bins with no presences
predBins <- predBins[predBins$Freq > 0,]

### Scale covariates

# Scale the prediction steps for broadleaf and coniferous woodland separately, and connectivity
predBins$coverScaled <- ( predBins$Cover - 
                               scalingParams[scalingParams$variable == "coverBF", "variableMean"] ) /
  scalingParams[scalingParams$variable == "coverBF", "variableSD"]

predBins$connScaled <- 
  (predBins$Connectivity - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
  scalingParams[scalingParams$variable == "connW", "variableSD"]

# Calculate scaled interaction terms for prediction
predBins$coverConnInt <- predBins$coverScaled * predBins$connScaled

# Calculate simplified prediction from coverScaled, connectivity, and their interaction
predBins$Prediction <- predBins$coverScaled  * 
  mean(spread_draws(coverWood_brms, b_Intercept)$b_Intercept)   +
  predBins$connScaled  *  
  mean(spread_draws(connWood_brms, b_Intercept)$b_Intercept) +
  predBins$coverConnInt * 
  mean(spread_draws(intWood_brms, b_Intercept)$b_Intercept)

# Plot
ggplot() +
  ggtitle("Woodland species") +
  xlab("Connectivity") +
  ylab("Relative occurence probability") +
  scale_colour_gradient(low = "red", high = "green", na.value = NA) +
  labs(colour =  "Cover") +
  theme_minimal() +
  geom_line(aes(x = predBins$Connectivity,
                y = 1 - exp(-exp(predBins$Prediction)), #1 - exp(-exp(predBins$Prediction)) is cloglog link
                group = predBins$Cover ,
                colour = predBins$Cover), size = 1)


