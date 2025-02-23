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

# Set supplementary  directories
SI_dirs <-  c(
  paste0(dataDir, "Supplementary_below30_analysis/"),
  paste0(dataDir, "Supplementary_quadratic_analysis/"),
  paste0(dataDir, "Supplementary_quartile_analysis/0-25_percent_cover_quartile/"),
  paste0(dataDir, "Supplementary_quartile_analysis/25-50_percent_cover_quartile/"),
  paste0(dataDir, "Supplementary_quartile_analysis/50-75_percent_cover_quartile/"),
  paste0(dataDir, "Supplementary_quartile_analysis/75-100_percent_cover_quartile/")
  )

# Set names of directories
names(SI_dirs) <- c("SI_below_30",
                    "SI_quadratic",
                    "SI_quartile_0-25",
                    "SI_quartile_25-50",
                    "SI_quartile_50-75",
                    "SI_quartile_75-100"
                    ) 

# Load main SDM fixed effect summaries
load(file = "../Data/Species_data/SDM_fixed_effect_summaries.RData")

# SET UP LOOP THROUGH SUPPLEMENTARY ANALYSES -----------

for(SI_dir in names(SI_dirs)) {

# SET UP LOOP THROUGH TAXA GROUPS ----------------------

# Create empty data frame to populate with all species level effects
allSpEff_df <- data.frame()

# Loop through each taxa group here
for (i in taxaGroups) {

  # Set taxa group
  taxaGroup <- i; print(taxaGroup)
  
  # Set taxa group folder
  taxaDir <- paste0(SI_dirs[SI_dir], taxaGroup)
  
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
    # (Remove file name and directory path)
    iSpeciesEffects$Species <-
      sub("/modelSummary.RData",
          "",
          allSummaryFiles[j]) %>%
      sub(paste0('.*', taxaGroup, '/'), 
          '',
          .)
      
    # Load species 'j' summary file into global environment
    grep(paste0(iSpeciesEffects$Species, "/modelSummary.RData"),
         fixed = TRUE,
         allSummaryFiles,
         value = TRUE) %>% 
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
  
  # Add to allSpEff_df dataframe
  allSpEff_df <- bind_rows(allSpEff_df, speciesEffects_df)
  
}

# CHECK DUPLCIATES -----------------------------------------------

# Find all duplicate rows
duplicates <- allSpEff_df %>% 
  group_by(species) %>% 
  filter(n()>
           ( NROW(speciesEffects_df) / length(allSummaryFiles) ) 
            )# 8 effect sizes per species

# Print if any duplicates - should be 0
print(duplicates)

# RESTRUCTURE DATA FRAME ------------------------------

# Drop unneeded columns, then spread dataframe so that each effect mean, sd, and
# quantile is a separate column
SI_df <-
  dplyr::select( allSpEff_df, -c("mode", "kld", "X0.5quant", )) %>%
  pivot_wider(
    names_from = effect,
    values_from = c( "mean", "sd", "X0.025quant", "X0.975quant" ))

# Remove any species with 0s for all mean and sd values (modelling error)
SI_df <- SI_df %>%
  filter(rowSums(abs(SI_df[-c(1:2)])) > 0)

# IDENTIFY WOODLAND AND CONNECTIVITY RELATIONSHIPS --------

# Use main analysis to assign relationships
# N.B. Inner join only retains species from main analysis for consistency
SI_df <- inner_join(
  meta_df %>%
    select(
      species,
      broadleafAssociation,
      coniferousAssociation,
      openAssociation,
      connectivitySig
      ), 
  SI_df,
  by = "species")

# SAVE DATA FRAME -----------------------------------------------

save(SI_df,
     file = paste0("../Data/Species_data/SDM_fixed_effect_summaries_",
                   SI_dir,
                   ".RData"))
}
