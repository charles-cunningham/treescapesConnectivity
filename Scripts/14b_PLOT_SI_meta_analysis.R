# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Plot meta analysis
#
# Script Description: Plot figures showing results of brms meta-analysis

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
library(ggridges)
library(rphylopic)
library(wesanderson)

### RELOAD OBJECTS ---------------------------------------------------

# List of brms objects
brmsList <- list(c("connBF_brms_under30", "coverBF_brms_under30", "intBF_brms_under30"),
              c("connBF_brms_quad", "coverBF_brms_quad", "quadBF_brms_quad", "intBF_brms_quad"),
              c("connBF_brms_0_25"), 
              c("connBF_brms_25_50"),
              c("connBF_brms_50_75"),
              c("connBF_brms_75_100"))

names(brmsList) <- c("under_30", "quadratic",
                     "quartile_0_25", "quartile_25_50", "quartile_50_75", "quartile_75_100")

# Load brms objects
load(file = "../Data/Species_data/brms_objects_below_30.RData")
load(file = "../Data/Species_data/brms_objects_quadratic.RData")
load(file = "../Data/Species_data/brms_objects_quart.RData")

# SET PARAMETERS ------------------------------------

# Set taxa groups to analyse
taxaGroups <- c( "Butterflies",  "Caddisflies", "Carabids", 
                 "Centipedes", "Ephemeroptera", "Gelechiidae", 
                 "Hoverflies", "Ladybirds", "Molluscs",
                 "Moths", "Odonata", "Orthoptera",
                 "Shieldbugs", "Soldierflies", "Spiders")

# Set taxa group labels
taxaGroupLabels <- c( "Butterflies", "Caddisflies", "Ground beetles",
                      "Centipedes", "Mayflies", "Gelechiid moths",
                      "Hoverflies","Ladybirds", "Molluscs",
                      "Moths", "Dragonflies", "Crickets/Grasshoppers",
                      "Shieldbugs", "Soldierflies", "Spiders")

# Set taxa group labels
taxaGroupLabels<- c( "Butterflies" = "Butterflies",
                     "Caddisflies" = "Caddisflies",
                     "Carabids" = "Ground beetles",
                     "Centipedes" = "Centipedes",
                     "Ephemeroptera" = "Mayflies",
                     "Gelechiidae" = "Gelechiid moths",
                     "Hoverflies" = "Hoverflies",
                     "Ladybirds" = "Ladybirds",
                     "Molluscs" = "Molluscs",
                     "Moths" = "Moths",
                     "Odonata" = "Dragonflies",
                     "Orthoptera" = "Crickets/Grasshoppers",
                     "Shieldbugs" = "Shield bugs",
                     "Soldierflies" = "Soldierflies",
                     "Spiders" = "Spiders")

# Set phylopic images (choose uuid manually)
phylopicImages <- data.frame(taxa = taxaGroups,
                             uuid = c("5aeaf558-3c48-4173-83b4-dbf2846f8d75",
                                      "d04ff776-292a-4181-a7d9-117a7036a348",
                                      "7832cac2-f122-4112-961e-2d506789e1c1",
                                      "fd2e5ed5-c6ae-4c84-ad86-db6620e7f967",
                                      "e642c1c1-76c1-4806-a024-aa9737d5bc41",
                                      "fae32bfc-d417-411d-be36-aea22c4a4a07",
                                      "fb360445-4c23-452e-b43d-34b42ce449dd",
                                      "984448b7-2ada-4e49-aaf1-3dcb7b056532",
                                      "cdc067a4-c2cc-4962-b330-2c9d02d245cb",
                                      "4afa48c4-0147-4b6c-b206-eb601c821a65",
                                      "2ef6d810-dc60-4dd2-aaba-bc0313c66ff6",
                                      "7c142ec5-aebb-495d-80fa-1b575090d5db",
                                      "bb330abc-f3ed-4d0e-a419-6fecf71ece1a",
                                      "8fa77803-7f94-4f0c-b2fd-66c2f70a03b4",
                                      "c27d24d2-7c1e-4dbd-8f0b-d2f7009bbe3c")) %>%
  mutate(svg = lapply(uuid, get_phylopic)) # Use uuids to get image objects

# Get attribution using get_attribution():
#get_attribution(uuid = "5aeaf558-3c48-4173-83b4-dbf2846f8d75")

# Get standardised species from quadratic analyses (the ones in the 75-100% quartile)
load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quartile_75-100.RData")
quadSpecies <- unique(SI_df$species)
rm(SI_df)

# SET UP LOOP THROUGH SI HERE -------------------------

# Loop through every SI analysis
for (SI in names(brmsList)) {

  # Load SDM fixed effect summaries from analysis ( and filter to quadSpecies if quartile)
  if (SI == "under_30") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_below_30.RData")
    
  } else if (SI == "quadratic") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quadratic.RData")
    
  } else if (SI == "quartile_0_25") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quartile_0-25.RData")
    SI_df <- SI_df %>%
      filter(species %in% quadSpecies)
    
  } else if (SI == "quartile_25_50") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quartile_25-50.RData")
    SI_df <- SI_df %>%
      filter(species %in% quadSpecies)
    
  } else if (SI == "quartile_50_75") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quartile_50-75.RData")
    SI_df <- SI_df %>%
      filter(species %in% quadSpecies)
    
  } else if (SI == "quartile_75_100") {
    
    load(file = "../Data/Species_data/SDM_fixed_effect_summaries_SI_quartile_75-100.RData")
  }
  
# SUMMARISE EFFECT SIZES -------------------------------

# Summarise effect sizes by  woodland association and taxa
summary_taxa_df <- SI_df %>%
  # Group by taxa group and effect category...
  group_by(connectivitySig, taxa,
           broadleafAssociation, coniferousAssociation, openAssociation) %>%
  # ... and count total number
  summarise(nuSpecies = length(species)) %>%
  # Ungroup for ...
  ungroup %>%
  # ... proportion
  mutate(freq = prop.table(nuSpecies),
         .by = taxa, connectivitySig) %>%
  # # ...and total species per taxa
  mutate(nuTaxa = sum(nuSpecies), .by = taxa)

# Summarise effect sizes by woodland association
summary_pooled_df <- SI_df %>%
  # Group by taxa group and effect category...
  group_by(connectivitySig,
           broadleafAssociation, coniferousAssociation, openAssociation ) %>%
  # ... and count total number
  summarise(nuSpecies = length(species)) %>%
  # Ungroup for ...
  ungroup %>%
  # ... proportion
  mutate(freq = prop.table(nuSpecies),
         .by = connectivitySig) %>%
  # Add taxa column
  add_column(taxa = "Pooled species") %>%
  # Add total number of species column
  mutate(nuTaxa = sum(nuSpecies))

# Join summary data frames together
summary_all_df <- rbind(summary_taxa_df,
                        summary_pooled_df)

# BASIC ANALYSIS ---------------------------------------------

### TAXA SUMMARIES

# Loop through brms object list
for (i in brmsList[SI][[1]]) {
  
  # Assign brms object
  iBrms <- get(i)
  
  ### PROCESS BAYESIAN DRAWS 
  
  # Extract draws grouping by taxa
  studyDraws <- spread_draws(iBrms,
                             r_taxa[taxa, ],
                             b_Intercept) %>%
    mutate(taxa_mean = r_taxa + b_Intercept)
  
  # Extract draws (all pooled)
  pooledEffectDraws <- spread_draws(iBrms, b_Intercept) %>%
    rename(taxa_mean = b_Intercept) %>%
    mutate(taxa = "Pooled species")
  
  # Bind taxa draws and pooled draws together, and reorder
  forestData <- bind_rows(studyDraws,
                          pooledEffectDraws) %>%
    
    ungroup() %>%
    mutate(taxa = reorder(taxa, taxa_mean)) %>%
    mutate(taxa = relevel(taxa,
                          "Pooled species",
                          after = Inf))
  
  # Group by taxa and summarize by intercept
  forestDataSummary <- group_by(forestData, taxa) %>%
    mean_qi(taxa_mean)
  
  # Print taxa summaries, and assign out for plotting
  print( i )
  print( forestDataSummary )
  assign(paste0(i, "_draws"), forestData)
  assign(paste0(i, "_drawsSummary"), forestDataSummary)
}

# PLOT POOLED ESTIMATES ------------------------------------------

# Loop through meta analysis subsets
for (i in brmsList[SI][[1]]) {
  
  # Assign brms objects
  iBrms <- get(i)
  iBrmsDraws <- get(paste0(i, "_draws"))
  iBrmsSummary <- get(paste0(i, "_drawsSummary"))
  
  # Set x label for taxa silhouettes for each brms run
  if (i %in% c(
    "connBF_brms_under30",
    "connBF_brms_quad",
    "connBF_brms_0_25",
    "connBF_brms_25_50",
    "connBF_brms_50_75",
    "connBF_brms_75_100"
  )) {
    xLabel <- "Connectivity effect estimate"
  } else if (i %in% c("coverBF_brms_under30", "coverBF_brms_quad")) {
    xLabel <- "Cover effect estimate"
  } else if (i %in% c("intBF_brms_under30", "intBF_brms_quad")) {
    xLabel <- "Cover-connectivity interaction effect estimate"
  } else if (i %in% c("quadBF_brms_quad")) {
    xLabel <- "Cover quadratic effect estimate"
  }
  
  ### PLOT 
  
  taxaSummaries <- ggplot(data = iBrmsDraws,
                          aes(taxa_mean,
                              taxa,
                              fill = taxa )) +
    
    # Add densities
    geom_density_ridges(scale = 0.95,
                        rel_min_height = 0.01) +
    geom_pointinterval(data = iBrmsSummary,
                       linewidth = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Change colours and labels
    scale_y_discrete(labels = taxaGroupLabels) +
    
    # Add vertical lines for pooled effect mean and CI, and 0
    geom_vline(xintercept = fixef(iBrms)[1, 1],
               color = "grey",
               linewidth = 1) +
    geom_vline( xintercept = fixef(iBrms)[1, 3:4],
                color = "grey",
                linetype = 2) +
    geom_vline(xintercept = 0,
               color = "black",
               linewidth = 1) +
    
    # Add taxon silouettes
    geom_phylopic(data = data.frame(taxa = levels(iBrmsDraws$taxa )) %>%
                    left_join(., phylopicImages , by = "taxa"),
                  inherit.aes = FALSE,
                  aes(x = min(iBrmsDraws$taxa_mean) +
                        (max(iBrmsDraws$taxa_mean) - 
                           min(iBrmsDraws$taxa_mean)) / 100,
                      y = taxa,
                      img = svg ),
                  size = 0.7,
                  na.rm = TRUE) +
    
    # Add labels
    labs(x = xLabel, # summary measure
         y = element_blank()) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16))
  
  # Save
  ggsave(filename = paste0("../Writing/Plots/SI/", "Meta_",
                           i, "_taxaSummaries.png"),
         taxaSummaries,
         dpi = 600,
         units = "px",
         width = 5000,
         height = 5000)
}

# PLOT POOLED ESTIMATES AND INDIVIDUAL SPECIES EFFECTS -------------------------

# Loop through SI analysis
for (i in brmsList[SI][[1]]) {

# Select out connectivity brms from within SI analysis
if (i %in% c(
  "connBF_brms_under30",
  "connBF_brms_quad",
  "connBF_brms_0_25",
  "connBF_brms_25_50",
  "connBF_brms_50_75",
  "connBF_brms_75_100"
)) {
  
  # Assign brms object
  iBrms <- get(i)
  iBrmsDraws <- get(paste0(i, "_draws"))
  iBrmsSummary <- get(paste0(i, "_drawsSummary"))
  
  ### PLOT 
  
  taxaSummaries <- ggplot(data = iBrmsDraws,
                          aes(taxa_mean,
                              taxa,
                              fill = taxa )) +
    
    # Add densities
    geom_density_ridges(scale = 0.95,
                        rel_min_height = 0.01) +
    geom_pointinterval(data = iBrmsSummary,
                       linewidth = 2,
                       aes(xmin = .lower, xmax = .upper)) +
    
    # Change colours and labels
    scale_y_discrete(labels = taxaGroupLabels) +
    
    # Add vertical lines for pooled effect mean and CI, and 0
    geom_vline(xintercept = fixef(iBrms)[1, 1],
               color = "grey",
               linewidth = 1) +
    geom_vline( xintercept = fixef(iBrms)[1, 3:4],
                color = "grey",
                linetype = 2) +
    geom_vline(xintercept = 0,
               color = "black",
               linewidth = 1) +
    
    # Add labels
    labs(x = "Connectivity effect estimate", # summary measure
         y = element_blank()) +
    
    # Add theme elements
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          plot.title = element_text(size = 22, hjust = -0.4, vjust =-0.2))
  
  ### BAR PLOT 
  
  # Set relevant data frame, and 
  barData <- subset(summary_all_df,
                    broadleafAssociation == "Y")
  
  # Set x position for taxa silhouettes
  if(i == "connBF_brms_under30") { 
    phyloX <- 450
  } else if (i == "connBF_brms_quad") { 
    phyloX <- 450
  } else if (i == "connBF_brms_0_25"){ 
    phyloX <- 200
  } else if (i == "connBF_brms_25_50"){ 
    phyloX <- 200
  } else if (i == "connBF_brms_50_75"){ 
    phyloX <- 200
  } else if (i == "connBF_brms_75_100"){ 
    phyloX <- 200
  }

  # Create bar plot
  sigBarPlot <- ggplot(data = barData,
                       aes(x = factor(taxa,
                                      levels = levels(iBrmsDraws$taxa)), 
                           y = nuSpecies,
                           fill = factor(connectivitySig,
                                         levels = c("Pos", "NS", "Neg")))) +
    
    # Add bar plot
    geom_bar(position = "stack",
             stat = "identity") +
    coord_flip() +
    
    # Change scales, fills and labels
    scale_x_discrete(labels = NULL) +
    scale_fill_manual(
      "Connectivity effect",
      values = wes_palette("Zissou1")[c(1,3,5)],
      labels = c("Positive", "None", "Negative"),
      guide_coloursteps(title.position = "top")) +
    labs(y = "Number of species") +
    
    # Add phylopic images
    geom_phylopic(data = data.frame(taxa = levels(iBrmsDraws$taxa)) %>%
                    left_join(., phylopicImages,
                              by = "taxa"),
                  inherit.aes = FALSE,
                  aes(x = taxa,
                      y = phyloX,
                      img = svg),
                  position = position_nudge(y = 0.5),
                  size = 0.8,
                  na.rm = TRUE) +
    
    # Change theme parameters
    theme_classic() +
    theme(axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          panel.grid.major.x= element_line( linewidth = 0.1, color = "black" ),
          legend.position = c(0.8, 0.8),
          legend.background = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          plot.title = element_text(size=22, hjust = -0.4, vjust =-0.2))
  
  ### COMBINE AND SAVE
  
  combinedPlot <- cowplot::ggdraw(clip = "on") +
    cowplot::draw_plot(taxaSummaries, 0, 0, 0.575, 1) +
    cowplot::draw_plot(sigBarPlot, 0.575, 0, 0.425, 1,) +
    cowplot::draw_label("(a)", 0.015, 0.98, size = 24) +
    cowplot::draw_label("(b)", 0.56, 0.98, size = 24) +
    theme(plot.background = element_rect( fill = "white", colour = "white"))
  
  # Save
  # Error message like "Removed X rows containing non-finite values (`stat_density_ridges()`)"
  # is expected due to manual x limit setting for broadleaf connectivity plot
  ggsave(filename = paste0("../Writing/Plots/SI/", "Meta_",
                           i, "_taxaMultiPlot.png"),
         combinedPlot,
         dpi = 600,
         units = "px",
         width = 10000,
         height = 5000)
}
}
}