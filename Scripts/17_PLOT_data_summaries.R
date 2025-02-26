# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: SDM spatial analysis
#
# Script Description: Analysis of SDM inlabru model outputs to investigate 
# potential spatial priorities for targeting improving connectivity

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(GGally)

# SET PARAMETERS ------------------------------------

# Set taxa group labels
taxaGroupLabels <- c( "Butterflies" = "Butterflies",
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

# DATA FILES ------------------------------------------

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# SpatRasters
for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                     i))) 
}

# DATA SUMMARIES --------------------------------------

# Proportion of converged models which are broadleaf- coniferous- ans open-associated species
NROW(subset(meta_df,broadleafAssociation == "Y")) / NROW(meta_df) * 100
NROW(subset(meta_df, coniferousAssociation == "Y")) / NROW(meta_df) * 100
NROW(subset(meta_df, openAssociation == "Y")) / NROW(meta_df) * 100

# Proportion of species which are both broadleaf and coniferous associated
NROW(subset(meta_df,broadleafAssociation == "Y" & coniferousAssociation == "Y"))
NROW(subset(meta_df,broadleafAssociation == "N" & coniferousAssociation == "N" & openAssociation == "N"))

# Proportion of broadleaf-associated species which have positive/negative/no connectivity effect
NROW(subset(meta_df,broadleafAssociation == "Y" & connectivitySig == "Pos")) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100
NROW(subset(meta_df,broadleafAssociation == "Y" & connectivitySig == "Neg")) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100
NROW(subset(meta_df,broadleafAssociation == "Y" & connectivitySig == "NS")) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100

# Proportion of species with positive/negative/no cover:connectivity interaction effect
NROW(subset(meta_df,broadleafAssociation == "Y" & X0.025quant_BFconnINT > 0 & X0.975quant_BFconnINT > 0)) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100
NROW(subset(meta_df,broadleafAssociation == "Y" & X0.025quant_BFconnINT < 0 & X0.975quant_BFconnINT < 0)) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100
NROW(subset(meta_df,broadleafAssociation == "Y" & X0.025quant_BFconnINT <= 0 & X0.975quant_BFconnINT >= 0)) / 
  NROW(subset(meta_df,broadleafAssociation == "Y")) * 100

# BROADLEAF COVER HISTOGRAM ------------------------

# Find quants
quants <- coverBF$BF_1990[] %>%
  .[!is.na(.)] %>%
  .[. > 0] %>%
  quantile(., probs = seq(0,1,0.25)) 

# Plot
coverHist <- ggplot(data.frame(coverBF)) +
  geom_histogram(aes(x = BF_1990), binwidth = 0.01) +
  geom_vline(xintercept = quants,
             color = "blue",
             size = 0.1) +
  theme_classic() +
  labs(x = "Proportion broadleaf cover", y = "Count")

#Save
ggsave(filename = paste0("../Writing/Plots/",
                         "coverHistogram.png"),
       coverHist,
       dpi = 600,
       units = "px",
       width = 3000,
       height = 2000)


# EFFECTS CORRELATION MATRIX ------------------------------------------

# Change taxa names to standardised names
# (as tricky to use labeller within ggpairs() as you normally would with ggplot())
meta_df$labels <- str_replace_all(string = meta_df$taxa,
                                  pattern = taxaGroupLabels)

# Plot
effectPairs <- 
  # Subset to broadleaf species
  subset(meta_df, broadleafAssociation == "Y") %>%
  # Select columns we want for pairs plot
  select(c(labels, mean_coverBF, mean_connectivity, mean_BFconnINT)) %>%
  # Call ggpairs
  ggpairs(
    aes(colour = labels, alpha = 0.5),
    lower = list(continuous = wrap("cor", size = 3)),
    diag = list(continuous = wrap("densityDiag")),
    upper = list(
      combo = wrap("box_no_facet"),
      continuous = wrap("points")
    ),
    columnLabels = c(
      "Recording scheme",
      "Broadleaf cover effect",
      "Connectivity effect",
      "Interaction effect"
    ))

# Remove 'Taxa' histograms as redundant information and confusing

# Create plots list
plots = list()

# Iteratively collect the sub-plots we want
for (i in 1:4) {
  plots <- c(plots, lapply(2:effectPairs$ncol, function(j)
    getPlot(
      effectPairs, i = i, j = j
    )))
}

# Use ggmatrix to plot the subsetted sub-plots
effectPairsTrim <- ggmatrix(
  plots,
  nrow = 4,
  ncol = effectPairs$ncol - 1,
  xAxisLabels = effectPairs$xAxisLabels[2:effectPairs$ncol],
  yAxisLabels = effectPairs$yAxisLabels
)

# Save
ggsave(
  filename = paste0("../Writing/Plots/", "rawPairsPlot.png"),
  effectPairsTrim,
  dpi = 600,
  units = "px",
  width = 6000,
  height = 6000
)
