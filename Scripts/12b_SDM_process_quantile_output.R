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

# Set taxa groups to analyse (only these three had converged models for all quantiles)
taxaGroups <- c( "Butterflies", "Moths", "Caddisflies")

# Set quantiles
quants <- c("0-20", "20-40", "40-60", "60-80", "80-100")

# DATA FILES ------------------------------------------

# SDM data folder
dataDir <- paste0("../Data/Species_data/Supplementary_cover_analysis/")

# Load SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# SET UP LOOP THROUGH TAXA GROUPS ----------------------

# Create empty data frame to populate with all species level effects
allSpEff_df <- data.frame()

# Loop through each taxa group here
for (i in taxaGroups) {
 
  # Set taxa group
  taxaGroup <- i; print(taxaGroup)
  
  # Set taxa group folder
  taxaDir <- paste0(dataDir, taxaGroup)
  
  # Find converged species (ones where last quantile [80-100] exists)
  covergedSpecies <- paste0(taxaDir, "/80-100_percent_cover_quantile") %>%
    list.files(.,
               recursive = TRUE) %>%
    gsub( "\\/.*","", . ) %>%
    unique()
  
  # Loop through each quantile here
  for (quant in quants) {

    # Set quantile directory
    quantDir <- paste0(taxaDir, "/", quant, "_percent_cover_quantile")
    
    # List summary files
    allSummaryFiles <- list.files(quantDir, 
                                  full.names = TRUE, 
                                  recursive = TRUE) %>%
      .[grep("modelSummary.RData", .)]
    
    # Subset to only converged species
    allSummaryFiles <- lapply(covergedSpecies, function(x){
      
      grep(x, allSummaryFiles)
      
      } ) %>% 
      unlist() %>%
      allSummaryFiles[.]
    
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
        sub(paste0(dataDir, taxaGroup, "/"),
            "",
            allSummaryFiles[j]) %>%
        sub("/modelSummary.RData", "", .)
      
      # Load species 'j' summary file into global environment
      grep(
        paste0(iSpeciesEffects$Species, "/modelSummary.RData"),
        fixed = TRUE,
        allSummaryFiles,
        value = TRUE
      ) %>%
        load(., envir = .GlobalEnv)
      
      # Assign fixed effects to list
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
    
    # Add quant column
    speciesEffects_df$quant <- quant
    
    # Add to allSpEff_df dataframe
    allSpEff_df <- bind_rows(allSpEff_df, speciesEffects_df)
    
  }}
  
  # CHECK DUPLCIATES -----------------------------------------------
  
  # Find all duplicate rows
  duplicates <- allSpEff_df %>%
    group_by(species, quant) %>%
    filter(n() > 8) # 8 effect sizes per species
  
  print(duplicates)
  
  # RESTRUCTURE DATA FRAME ------------------------------
  
  # Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
  # quantile is a separate column
  quant_df <-
    dplyr::select(allSpEff_df, -c("mode", "kld", "X0.5quant", )) %>%
    pivot_wider(
      names_from = effect,
      values_from = c("mean", "sd", "X0.025quant", "X0.975quant")
    )
  
  # Remove any species with 0s for all mean and sd values (modelling error)
  quant_df <- quant_df %>%
    filter(rowSums(abs(quant_df[c(4:19)])) > 0)
  
  # SAVE DATA FRAME -----------------------------------------------
  
  save(quant_df, file = "../Data/Species_data/Quant_fixed_effect_summaries.RData")
  