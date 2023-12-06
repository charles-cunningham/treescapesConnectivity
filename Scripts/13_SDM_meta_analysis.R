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
group_df <- meta_df %>%
  group_by( woodSp, connectivity_sig ) 

# Get basic numbers
summarise(group_df, length(species))

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

# GROUPED COVER AND CONNECTIVITY ANALYSIS -----------------------------------

# Create raster for sum of every pooled group

# Loop through every group in group_df
groupSumsR <- lapply(1:NROW(group_keys(group_df)) , function(i) {
  
  # Split group into separate df, get df i...
  iGroup <- group_split(group_df)[[i]] %>%
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
  names(iGroupSum) <- group_keys(group_df)[i, ] %>%
    as.character  %>%
    paste(., collapse = " ")
  
  return(iGroupSum)
  
# Join all group sum spatRasts into single spatRast  
}) %>% rast(.) 

# Calculate mean
groupMeansR <- groupSumsR /
  sapply(group_split(group_df), NROW)

# Save
writeRaster(groupSumsR,
            paste0("../Data/Species_data/Summed_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)
writeRaster(groupMeansR,
            paste0("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif"),
            overwrite = TRUE)

### Plot

# Create plot directory
dir.create("../Writing/Plots/GroupedOccPlots")

# Loop through each cover/connectivity group
for(i in 1:nlyr(groupMeansR)) {

  # Create data frame from grouped mean spatRast i
  iGroupMeans_df <- groupMeansR[[i]] %>%
    as.data.frame(., xy = TRUE)
  
  # GGplot
  occPlot <- ggplot(data = iGroupMeans_df,
                    aes(x = x, y = y,
                      colour = iGroupMeans_df[, 3],
                      fill = iGroupMeans_df[, 3])) +
    
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
      label = paste(NROW(group_split(group_df)[[i]]),
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
                           names(iGroupMeans_df[3]), ".png"),
    occPlot,
    dpi = 600,
    units = "px",
    width = 6000,
    height = 8000)

}

# INTERACTION PLOTS ----------------------------------- IN PROGRESS, temp version for now

load("../Data/Spatial_data/DataForInlabru/scalingParams.RData")

# How many prediction steps?
nSamp <- 5

# Extract max scaled value
maxConnectivity <- rast("../Data/Spatial_data/DataForInlabru/spatRaster/connW.tif") %>%
  global(., fun = "max", na.rm = TRUE) %>% 
  max

### Set up prediction data frames

# Create unscaled data frame of cover and connectivity values to predict over
# Seperate broadleaf and coniferous data frames
BF_pred_df <-  expand.grid(BF_pred = seq(0, 1, by = 1/nSamp),
                           conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))

### Scale covariates

# Scale the prediction steps for broadleaf and coniferous woodland seperately, and connectivity
# N.B. Have to name columns the same as the original SGDF datasets!
BF_pred_df$coverBF_SGDF <- ( BF_pred_df$BF_pred - 
                               scalingParams[scalingParams$variable == "coverBF", "variableMean"] ) /
  scalingParams[scalingParams$variable == "coverBF", "variableSD"]

BF_pred_df$connW_SGDF <- 
  (BF_pred_df$conn_pred - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
  scalingParams[scalingParams$variable == "connW", "variableSD"]

# Calculate scaled interaction terms for prediction
BF_pred_df$coverBF_connW_SGDF <- BF_pred_df$coverBF_SGDF * BF_pred_df$connW_SGDF


group_keys(group_df)
test <- group_split(group_df)[[9]]
test <- rbind(group_split(group_df)[[7]], group_split(group_df)[[8]], group_split(group_df)[[9]])
#test <- group_df

BF_pred_df$Prediction <- BF_pred_df$coverBF_SGDF  * mean(test$mean_coverBF)    +
  BF_pred_df$connW_SGDF  *  mean(test$mean_connectivity ) +
  BF_pred_df$coverBF_connW_SGDF *  mean(test$mean_BFconnINT ) + 
  mean(test$mean_Intercept) 
  

ggplot() +
  ggtitle("Species with negative woodland association") +
  xlab("Connectivity") +
  ylab("Relative occurence probability (link scale)") +
  labs(colour =  "Cover") +
  theme_minimal() +
  geom_line(aes(x = BF_pred_df$conn_pred, y = BF_pred_df$Prediction, group = BF_pred_df$BF_pred ,
                colour = BF_pred_df$BF_pred))

#1 - exp(-exp(BF_pred_df$Prediction))

# CREATION PRIORITY MAPS ----------------------------- IN PROGRESS

groupSumsR <- rast("../Data/Species_data/Summed_occupancy_grouped_by_effects.tif")
groupMeansR <- rast("../Data/Species_data/Mean_occupancy_grouped_by_effects.tif")

test <- c(groupMeansR[[7]], groupMeansR[[8]], groupMeansR[[9]]) %>%
            mean %>%
  as.data.frame(., xy = TRUE)


ggplot(data = test) +
 geom_tile(aes(x=x,y=y, fill = mean, colour = mean)) +
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
  theme_void() + 
  coord_fixed() 








# Process BF cover
woodCover <- rast("../Data/Spatial_data/DataForInlabru/spatRaster/coverBF_scaled.tif")[["BF_2015"]]
woodCover <- project (woodCover, gsub( "units=km", "units=m",
                                       sf::st_crs(groupSumsR)$proj4string ))
crs(woodCover) <- crs(groupMeansR)

# Process connectivity
woodConn <- rast("../Data/Spatial_data/DataForInlabru/spatRaster/connW.tif")[["conn_2015"]]
woodConn <- project (woodConn, gsub( "units=km", "units=m",
                                       sf::st_crs(groupSumsR)$proj4string ))
crs(woodConn) <- crs(groupMeansR)



# Create bivariate data frame
bivariate_df <- c(woodConn, groupMeansR[["Prefer woodland Pos"]]) %>%
  as.data.frame(., xy = TRUE)

# Change names (no spaces, intuitive)
names(bivariate_df) <- c("x", "y", "Conn", "MeanOcc")

# Create 5 quantile buckets for occupancy
quantilesOcc <- bivariate_df %>%
  pull(MeanOcc) %>%
  quantile(probs = seq(0, 1, length.out = 5), na.rm = TRUE)
quantilesConn <- bivariate_df %>%
  pull(Conn) %>%
  quantile(probs = seq(0, 1, length.out = 5), na.rm = TRUE)

# Create 5x5 colour scale
bivariate_color_scale <- tibble(
  "1-1" = "#d3d3d3", # low x, low y
  "2-1" = "#b6cdcd",
  "3-1" = "#97c5c5",
  "4-1" = "#75bebe",
  "5-1" = "#52b6b6", # high x, low y
  "1-2" = "#cab6c5",
  "2-2" = "#aeb0bf",
  "3-2" = "#91aab9",
  "4-2" = "#70a4b2",
  "5-2" = "#4e9daa",
  "1-3" = "#c098b9",
  "2-3" = "#a593b3",
  "3-3" = "#898ead",
  "4-3" = "#6b89a6",
  "5-3" = "#4a839f",
  "1-4" = "#b77aab",
  "2-4" = "#9e76a6",
  "3-4" = "#8372a0",
  "4-4" = "#666e9a",
  "5-4" = "#476993",
  "1-5" = "#ad5b9c", # low x, high y
  "2-5" = "#955898",
  "3-5" = "#7c5592",
  "4-5" = "#60528d",
  "5-5" = "#434e87" # high x, high y
)%>%
 gather("group", "fill")

# Separate the groups
bivariate_color_scale <- bivariate_color_scale %>%
  separate(group, into = c("MeanOcc", "Conn"),
           sep = "-",
           remove = FALSE) %>%
  mutate(MeanOcc = as.integer(MeanOcc),
         Cover = as.integer(Conn))

# Cut into groups defined above and join fill
bivariate_df <- mutate(bivariate_df,
    Occ_quantiles = cut(MeanOcc,
                        breaks = quantilesOcc, # based on quantiles
                        include.lowest = TRUE),
    Conn_quantiles = cut(Conn,
                          breaks = 5, # equal interval cuts 
                          include.lowest = TRUE,
                          na.rm = TRUE),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste0(as.numeric(Conn_quantiles), "-",
                  as.numeric(Occ_quantiles))
  ) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")

### Create plot

# Create map
map <- ggplot(data = bivariate_df) +
  ggtitle(paste0("Modelled relative occurence proability of woodland species\n",
                 "with positive association with connectivity against woodland connectivity")) +
  geom_tile(aes(x=x,y=y, fill = fill, colour = fill),
            show.legend = FALSE) +
  scale_fill_identity() +
  scale_colour_identity() +
  theme_void() + 
  coord_fixed() 
  #geom_sf(data = st_as_sf(smoothUK_sp), fill = "NA", colour = "black", inherit.aes = FALSE)

# Create legend
legend <- ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = MeanOcc,
      y = Conn,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = "Woodland connectivity ⟶️",
       y = "Mean Occurence ⟶️") +
  theme_void() +
  # make font small enough
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, angle = 90) ) +
  # quadratic tiles
  coord_fixed()

# Plot map and legend
cowplot::ggdraw() +
  cowplot::draw_plot(map, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.65, 0.55, 0.25, 0.25)



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

