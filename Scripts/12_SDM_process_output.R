# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Extract SDM data
#
# Script Description: Pull out SDM inlabru model outputs

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(INLA)
library(inlabru)
library(tidyverse)
library(terra)

# SET PARAMETERS ------------------------------------

# Set taxa groups to analyse
taxaGroups <- c( "Butterflies", "Moths", "Carabids", 
              "Caddisflies", "Centipedes", "Ephemeroptera", 
              "Gelechiidae", "Hoverflies", "Ladybirds",
              "Molluscs", "Odonata", "Orthoptera",
              "Shieldbugs", "Soldierflies", "Spiders")

# DATA FILES ------------------------------------------

# SDM data folder
dataDir <- paste0("../Data/Species_data/SDMs/")

# SET UP LOOP THROUGH TAXA GROUPS ----------------------

# Create empty data frame to populate with all species level effects
allSpEff_df <- data.frame()

# Loop through each taxa group here
for (i in taxaGroups) {
  
  # Set taxa group
  taxaGroup <- i; print(taxaGroup)
  
  # Set taxa group folder
  taxaDir <- paste0(dataDir, taxaGroup)
  
  # List summary files
  allSummaryFiles <- list.files(taxaDir,
                                full.names = TRUE,
                                recursive = TRUE) %>%
    .[grep("modelSummary.RData",.)]

# IMPORT FIXED EFFECTS ---------------------------------

# Create empty list to populate with model fixed effect summaries
speciesEffects <- list()

# For every converged species model in taxa group...
for (j in 1:NROW(allSummaryFiles)) {

  # Create list just for that species (will become sub-list)
  iSpeciesEffects <- list()
  
  # Assign species name from file path to list
  # (Remove directory path and file name sequentially)
  iSpeciesEffects$Species <-
    sub(paste0("../Data/Species_data/SDMs/", taxaGroup, "/"),
        "",
        allSummaryFiles[j]) %>%
    sub("/modelSummary.RData",
        "",
        .)
  
  # Load species 'j' summary file into global environment
  grep(paste0(iSpeciesEffects$Species, "/modelSummary.RData"),
       fixed = TRUE,
       allSummaryFiles,
       value = TRUE) %>% 
    load(., envir = .GlobalEnv)
  
  # Assign fised effects to list
  iSpeciesEffects$FixedEffects <- data.frame(modelSummary$inla$fixed)
  
  # Save list of species i as a sub-list of speciesEffects
  speciesEffects[[j]] <- iSpeciesEffects
    
}

# Add in species name to fixed effects data frames,
# and effect from row names
for (j in 1:length(speciesEffects)) {
  speciesEffects[[j]][[2]]$species <- speciesEffects[[j]][[1]]
  speciesEffects[[j]][[2]]$effect <- rownames(speciesEffects[[j]][[2]])
}

# Extract fixed effects data frames and bind together
speciesEffects_df <- map(speciesEffects, 2, .default = NA) %>% 
  bind_rows

# Drop redundant row names
rownames(speciesEffects_df) <- NULL

# Add taxa column
speciesEffects_df$taxa <- taxaGroup

# Add to allSpEff_df dataframe
allSpEff_df <- bind_rows(allSpEff_df, speciesEffects_df)

}

# REMOVE DUPLCIATES -----------------------------------------------

# Find all duplicate rows
duplicates <- allSpEff_df %>% 
  group_by(species) %>% 
  filter(n()>8) # 8 effect sizes per species

# Loop to remove duplicate models
for (i in unique(duplicates$species)) {# For each duplicate species

  # Select models to drop based on output (drop species with most uncertainty)
  dropSpeciesRows <- duplicates %>% 
    filter(species == i ) %>% # Filter to species
    filter(effect == "connectivity" ) %>% # Filter to connectivity effect
    filter(sd != min(sd )) # Drop all but minimum sd 
  
  # Identify row numbers of models (connectivity effect only) to drop
  dropSpeciesRowN <- which(do.call(paste0, duplicates) %in%
                             do.call(paste0, dropSpeciesRows))
  
  # Find all rows to drop 
  # For each connectivity row number, get all effect rows (connectivity is the 4th row of 8)
  dropSpeciesAll <- lapply(dropSpeciesRowN, 
                          function(x) {(x-3):(x+4)}) %>% 
    unlist %>% # Merge into single vector of row numbers
    duplicates[.,] # Get associated rows

  # Drop rows from main data frame
  allSpEff_df <- anti_join(allSpEff_df, dropSpeciesAll)
}

# RESTRUCTURE DATA FRAME ------------------------------

# Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
# quantile is a separate column
meta_df <-
  dplyr::select( allSpEff_df, -c("mode", "kld", "X0.5quant", )) %>%
  pivot_wider(
    names_from = effect,
    values_from = c( "mean", "sd", "X0.025quant", "X0.975quant" ))

# IDENTIFY WOODLAND AND CONNECTIVITY RELATIONSHIPS --------

# Create woodland association categories for broadleaf and coniferous woodland:
# Four options: prefer broadleaf woodland only, prefer coniferous woodland only,
# prefer both woodland, prefer no woodland
meta_df <- mutate(meta_df,
                  broadleafAssociation = case_when(X0.025quant_coverBF   > 0 &
                                                     X0.975quant_coverBF > 0 ~ "Y",
                                                   TRUE ~ "N"),
                  coniferousAssociation = case_when(X0.025quant_coverCF   > 0 &
                                                      X0.975quant_coverCF > 0 ~ "Y",
                                                    TRUE ~ "N"),
                  openAssociation = case_when((X0.025quant_coverBF < 0 & X0.975quant_coverBF < 0 ) &
                                                (X0.025quant_coverCF < 0 & X0.975quant_coverCF < 0) ~ "Y",
                                              TRUE ~ "N"))

# Create discrete woodland association categories for connectivity
meta_df <- mutate( meta_df,
                   connectivitySig = case_when(
                     X0.025quant_connectivity > 0 & X0.975quant_connectivity > 0 ~ "Pos",
                     X0.025quant_connectivity < 0 & X0.975quant_connectivity < 0 ~ "Neg",
                     TRUE ~ "NS"))

# SAVE DATA FRAME -----------------------------------------------

save(meta_df,
     file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")
