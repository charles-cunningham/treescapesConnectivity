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
library(tidyverse)
library(tidybayes)
library(ggridges)

# SET PARAMETERS ------------------------------------

# Set taxa groups to analyse, and mobility class
taxa <- c( "Butterflies",
           "Moths",
           "Carabids",
           "Caddisflies",
           "Centipedes",
           "Ephemeroptera",
           "Gelechiidae",
           "Hoverflies",
           "Ladybirds",
           "Molluscs",
           "Odonata",
           "Orthoptera",
           "Shieldbugs",
           "Soldierflies",
           "Spiders" )

# Set placeholder mobility for each taxa group
mobility <- c( "high",
               "high",
               "medium",
               "high",
               "low",
               "medium",
               "high",
               "high",
               "high",
               "low",
               "high",
               "medium",
               "medium",
               "high",
               "medium" ) 

# Create data frame
taxaGroups_df <- data.frame(taxa, mobility)

# DATA FILES ------------------------------------------

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# RESTRUCTURE DATA FRAME ------------------------------

# Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
# quantile is a separate column
meta_df <-
  dplyr::select( allSpEff_df, -c("mode", "kld", "X0.5quant", )) %>%
  pivot_wider(
    names_from = effect,
    values_from = c( "mean", "sd", "X0.025quant", "X0.975quant" ))

# Add in mobility
meta_df <- inner_join(meta_df, taxaGroups_df,
                      by = join_by(taxa))

# IDENTIFY WOODLAND AND CONNECTIVITY REALATIONSHIPS -------------------

# Create discrete woodland association categories for BF and CF woodland
meta_df <- mutate( meta_df,
                   BFwoodSp = case_when(
                     X0.025quant_coverBF > 0 & X0.975quant_coverBF > 0 ~ "Prefer woodland",
                     X0.025quant_coverBF < 0 & X0.975quant_coverBF < 0 ~ "Avoid woodland",
                     TRUE ~ "Generalist"))
meta_df <- mutate( meta_df,
                   CFwoodSp = case_when(
                     X0.025quant_coverCF > 0 & X0.975quant_coverCF > 0 ~ "Prefer woodland",
                     X0.025quant_coverCF < 0 & X0.975quant_coverCF < 0 ~ "Avoid woodland",
                     TRUE ~ "Generalist"))
meta_df <- mutate( meta_df,
                   woodSp = case_when(
                     BFwoodSp == "Prefer woodland" | CFwoodSp == "Prefer woodland" ~ "Prefer woodland",
                     BFwoodSp == "Avoid woodland"  | CFwoodSp == "Avoid woodland"  ~ "Avoid woodland",
                     TRUE ~ "Generalist"))

# Create discrete woodland association categories for connectivity
meta_df <- mutate( meta_df,
                   connectivity_sig = case_when(
                     X0.025quant_connectivity > 0 & X0.975quant_connectivity > 0 ~ "Pos",
                     X0.025quant_connectivity < 0 & X0.975quant_connectivity < 0 ~ "Neg",
                     TRUE ~ "NS"))

# DESCRIPTIVE POOLED STATS

# Group by taxa group and effect category
grouped_df <- meta_df %>%
  group_by( woodSp, connectivity_sig ) 

# Get basic numbers
summarise(grouped_df, length(species))

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
allSummary <- meta_df %>%
  # Group by taxa group and effect category...
  group_by(taxa, connectivity_sig, woodSp ) %>% # Change to CFwoodSp to see coniferous associations
  # ... and count total number
  summarise(nuSpecies = length(species)) %>%
  # Ungroup for ...
  ungroup %>%
  # ... proportion
  mutate(freq = prop.table(nuSpecies), .by = c(taxa, woodSp))

# Plot
ggplot(allSummary) +
  geom_bar(aes(x = taxa, y = nuSpecies, fill = connectivity_sig),
           position = "stack",
           stat = "identity") +
  facet_wrap( ~ woodSp) + # Change to CFwoodSpp to see coniferous associations
coord_flip() 

ggplot(allSummary) +
  geom_bar(aes(x = taxa, y = freq, fill = connectivity_sig),
           position = "stack",
           stat = "identity") +
  facet_wrap( ~ woodSp) + # Change to CFwoodSpp to see coniferous associations
  coord_flip() 

# GROUPED COVER AND CONNECTIVITY ANALYSIS -----------------------------------

# Create raster for sum of every pooled group

# Loop through every group in grouped_df
groupSumsR <- lapply(1:NROW(group_keys(grouped_df)) , function(i) {
  
  # Split group into separate df, get df i...
  iGroup <- group_split(grouped_df)[[i]] %>%
    dplyr::select(species, taxa) #... and select taxa and species
  
  # Create list of file names of occupancy plots (Script 11) 
  # for every row in iGroup
  iGroupFiles <- apply(iGroup, 1, function(j) {
    paste0("../Data/Species_data/SDMs/",
           j["taxa"],
           "/",
           j["species"],
           "/occPlot.tif")
  })
  
  # Create spatRast
  iGroupRast <- rast(iGroupFiles)
  
  # Sum
  iGroupSum <- sum(iGroupRast)
  
  # Change spatRast layer name to group key
  names(iGroupSum) <- group_keys(grouped_df)[i, ] %>%
    as.character  %>%
    paste(., collapse = " ")
  
  return(iGroupSum)
  
# Join all group sum spatRasts into single spatRast  
}) %>% rast(.) 

# Calculate mean
groupMeansR <- groupSumsR /
  NROW(group_split(grouped_df))


# Save
writeRaster(groupSumsR,
            paste0("../Data/Species_data/Summmed_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)
writeRaster(groupMeansR,
            paste0("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)

### Plot

# Create plot directory
dir.create("../Writing/Plots/GroupedOccPlots")

# Loop through each cover/connectivity group
for(i in 1:nlyr(groupMeansR)) {
#for(i in 1:1) {


  # Create data frame from grouped sum spatRast i
  iGroupSum_df <- groupMeansR[[i]] %>%
    as.data.frame(., xy = TRUE)
  
  # GGplot
  occPlot <- ggplot(data = iGroupSum_df,
                    aes(x = x, y = y,
                      colour = iGroupSum_df[, 3],
                      fill = iGroupSum_df[, 3])) +
    
    # Set equal coordinates
    coord_equal() +
    
    # Country boundary polygons
    # geom_sf(data = st_as_sf(UK), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    # geom_sf(data = st_as_sf(Ireland), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    # geom_sf(data = st_as_sf(IsleOfMan), fill = "grey90", colour = NA, inherit.aes = FALSE) +
    
    # Rraster
    geom_tile() +
    scale_colour_gradient(low = "#e5f5e0",
                          high = "#31a354",
                          # limits = c(0,0.5),
                          guide = NULL) +
    scale_fill_gradient(
      "Mean relative\noccupancy probability",
      low = "#e5f5e0",
      high = "#31a354",
      #limits = c(0,0.5),
      guide = guide_colourbar(
        ticks = TRUE,
        draw.ulim = FALSE,
        draw.llim = FALSE,
        title.position = "top",
        label.position = "right",
        label.theme = element_text(size = 12),
        title.theme = element_text(size = 12),
        barwidth = unit(2, "lines"),
        barheight = unit(8, "lines") )) +
    
    # Add title
    annotate(
      geom = "text",
      x = 490000,
      y = 900000,
      label = names(iGroupMeans_df[3]),
      size = 6) +
    
    # Add title
    annotate(
      geom = "text",
      x = 90000,
      y = 300000,
      label = paste(NROW(group_split(grouped_df)[[i]]),
                    "species"),
      size = 6) +
    
    # Set theme parameters
    theme_void() +
    theme(plot.background = element_rect(fill = "white",
                                     colour = "white"),
          legend.position = c(0.8, 0.7),
          plot.margin = margin(-3, 0, -3, 0, "lines"))
  
  ggsave(filename = paste0("../Writing/Plots/GroupedOccPlots/",
                           "Occupancy_",
                           names(iGroupSum_df[3]), ".png"),
    occPlot,
    dpi = 600,
    units = "px",
    width = 6000,
    height = 8000)

}



# META ANALYSIS --------------------------------------
# In this stage we want to check overall data patterns,
# partitioning of variance between taxa and species levels etc.

### Connectivity

# Fit Bayesian meta analysis model
# N.B. Species effect nested within taxa effect.
# (i.e. taxa + taxa:species random effects)

# All species
metaAnAll <- brm(data = meta_df, 
              family = gaussian,
              mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 1), class = sd)),
              iter = 4000,
              cores = 4,
              chains = 4)

# Woodland species only
metaAnWood <- brm(data = meta_df[meta_df$BFwoodSp == "Prefer woodland" |
                                   meta_df$CFwoodSp == "Prefer woodland", ],
                 family = gaussian,
                 mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 1), class = sd)),
                 iter = 4000,
                 cores = 4,
                 chains = 4)

# Woodland species only
metaAnNotWood <- brm(data = meta_df[meta_df$BFwoodSp != "Prefer woodland" &
                                   meta_df$CFwoodSp != "Prefer woodland", ],
                  family = gaussian,
                  mean_connectivity | se(sd_connectivity) ~ 1 + (1 | taxa / species),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(cauchy(0, 1), class = sd)),
                  iter = 4000,
                  cores = 4,
                  chains = 4)

# List the different meta-analyses
metaAnList <- list(metaAnAll, metaAnWood, metaAnNotWood)

# Print summary for all species and woodland species
summary(metaAnAll)
summary(metaAnWood)
summary(metaAnNotWood)

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
  metaAn <- metaAnList[[i]]

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

### Cover:connectivity interaction

# Fit Bayesian meta analysis model with interaction
# N.B. Species effect nested within taxa effect.
# (i.e. taxa + taxa:species random effects)

# All species
metaAnIntAll <- brm(data = meta_df, 
                 family = gaussian,
                 mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 1), class = sd)),
                 iter = 4000,
                 cores = 4,
                 chains = 4)

# Woodland species only
metaAnIntWood <- brm(data = meta_df[meta_df$BFwoodSp == "Prefer woodland" |
                                   meta_df$CFwoodSp == "Prefer woodland", ],
                  family = gaussian,
                  mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(cauchy(0, 1), class = sd)),
                  iter = 4000,
                  cores = 4,
                  chains = 4)

# Woodland species only
metaAnIntNotWood <- brm(data = meta_df[meta_df$BFwoodSp != "Prefer woodland" &
                                      meta_df$CFwoodSp != "Prefer woodland", ],
                     family = gaussian,
                     mean_BFconnINT | se(sd_BFconnINT) ~ 1 + (1 | taxa / species),
                     prior = c(prior(normal(0, 1), class = Intercept),
                               prior(cauchy(0, 1), class = sd)),
                     iter = 4000,
                     cores = 4,
                     chains = 4)

# List the different meta-analyses
metaAnIntList <- list(metaAnIntAll, metaAnIntWood, metaAnIntNotWood)

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


##############################################################
# 
# ### TESTING
# 
# # META REGRESSION ------------------------------------
# # In this stage we check importance of different hypothesized 
# # predictors of connectivity importance
# 
# ### Cover (greater association with BL cover proxy for specialism)
# 
# # Initial plot
# meta_df %>%
#   filter(BFcoverMean < 10,
#          mean < 10,
#          mean > -20) %>% # lots of outliers
#   ggplot() +
#   geom_point(aes(x = BFcoverMean, y = mean, colour=taxa)) 
# 
# 
# # Fit Bayesian meta regression model inc. broadleaf cover
# metaRegCover <- brm(#data = meta_df,
#                     data = meta_df[meta_df$BFwoodSp == "Prefer woodland" |
#                                      meta_df$CFwoodSp == "Prefer woodland", ],
#               family = gaussian,
#               mean | se(sd) ~ 1 + BFcoverMean + (1 | taxa / species),
#               prior = c(prior(normal(0, 1), class = Intercept),
#                         prior(cauchy(0, 1), class = sd)),
#               iter = 4000,
#               cores = 4,
#               chains = 4)
# 
# 
# metaRegCover <- brm(data = meta_df[meta_df$BFwoodSp == "Prefer woodland", ],
#   family = gaussian,
#   mean | se(sd) ~ 1 + BFwoodSp + (1 | taxa / species),
#   prior = c(prior(normal(0, 1), class = Intercept),
#             prior(cauchy(0, 1), class = sd)),
#   iter = 4000,
#   cores = 4,
#   chains = 4)
# 
# 
# # Print summary
# summary(metaRegCover)
# 
# # Plot cover effect
# as_draws_df(metaRegCover) %>% 
#   
#   ggplot(aes(x = b_BFcoverMean, y = 0)) +
#   stat_halfeye(point_interval = median_qi, .width = c(.5, .95)) +
#   scale_y_continuous(NULL, breaks = NULL) +
#   theme(panel.grid = element_blank())
# 
# ### Mobility (taxa groups classified as 3 mobility classes)
# 
# # Initial plot
# meta_df %>%
#   filter(abs(mean) < 2.5) %>% # lots of outliers
#   ggplot() +
#   geom_boxplot(aes(x = mobility, y = mean))+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   coord_flip()
# 
# # Fit Bayesian meta regression model inc. mobility (factor)
# metaRegMob <- brm(#data = meta_df,
#                   data = meta_df[meta_df$BFwoodSp == "Prefer woodland" |
#                                    meta_df$CFwoodSp == "Prefer woodland", ],
#                family = gaussian,
#                mean | se(sd) ~ 1 + mobility + (1 | taxa / species),
#                prior = c(prior(normal(0, 1), class = Intercept),
#                          prior(cauchy(0, 1), class = sd)),
#                iter = 4000,
#                cores = 4,
#                chains = 4)
# 
# # Print summary
# summary(metaRegMob)
# 
# # Plot mobility effect
# as_draws_df(metaRegMob) %>% 
#   
#   ggplot(aes(x = b_mobilitylow, y = 0)) +
#   stat_halfeye(point_interval = median_qi, .width = c(.5, .95)) +
#   scale_y_continuous(NULL, breaks = NULL) +
#   theme(panel.grid = element_blank())
# 
# ### Mobility * specialism
# 
# # Fit Bayesian meta regression model inc. mobility (factor)
# metaRegInt <- brm(#data = meta_df,
#                   data = meta_df[meta_df$BFwoodSp == "Prefer woodland" |
#                                    meta_df$CFwoodSp == "Prefer woodland", ],
#                   family = gaussian,
#                   mean | se(sd) ~ 1 + coverEff * mobility + (1 | taxa / species),
#                   prior = c(prior(normal(0, 1), class = Intercept),
#                             prior(cauchy(0, 1), class = sd)),
#                   iter = 4000,
#                   cores = 4,
#                   chains = 4)
# 
# # Print summary
# summary(metaRegInt)
# 
# #####
# # testing
# ggplot(meta_df, aes(x=mean_connectivity)) + 
#   geom_histogram(color="black", fill="white", bins = 100) +
#   xlim(-5,5)
# 
# z <- meta_df[,c("mean_connectivity", "mean_coverBF", "mean_BFconnINT")] %>%
#   filter_all(all_vars(. < 5)) %>%
#   filter_all(all_vars(. > -5))
# 
# plot(z$mean_connectivity~ z$mean_coverBF)
# with(z,text(z$mean_connectivity~ z$mean_coverBF, labels=z$species,pos=1))
# 
# 
# means <- apply(z,2,mean)
# sds <- apply(z,2,sd)
# nor <- scale(z,center=means,scale=sds)
# 
# library(cluster)
# 
# kc<-kmeans(nor,5)
# kc$cluster
# z
# z <- cbind(z, cluster = kc$cluster)
# factoextra::fviz_cluster(kc, data = nor)
# 
# factoextra::fviz_nbclust(nor, kmeans, method = "wss")
# 
# ggplot(z) +
#   geom_point(mapping = aes(x = mean_connectivity, y = mean_coverBF, colour = as.character(cluster)))


