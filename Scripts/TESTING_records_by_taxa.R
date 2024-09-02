# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Species record comparison by taxa
#
# Script Description: Count total records for each taxonomic group for reporting
#

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)

# SOURCE FUNCTIONS ----------------------------

source("Functions/assignYearGroup().R")

# ASSIGN PARAMETERS ----------------------------------

# Organise raw data taxa groups
dataBC <- c( "Butterflies", "Moths" )
dataBRC <- c( "Bryophytes", "Carabids", "Caddisflies", "Centipedes",
              "Ephemeroptera", "Gelechiidae", "Hoverflies", "Ladybirds",
              "Lichen", "Molluscs", "Odonata", "Orthoptera",
              "Shieldbugs", "Soldierflies", "Spiders")

# List range of years used for species records
range1 = c(1990,2000)
range2 = c(2015,2025)

# ASSIGN TAXA RAW DATA FILE PATHS --------------------

# Create data frame of taxonomic groups and file path
taxaFiles <- 
  data.frame( 
    "taxaGroups" = c( "Butterflies", "Caddisflies", "Carabids",
                      "Centipedes", "Ephemeroptera", "Gelechiidae",
                      "Hoverflies", "Ladybirds","Molluscs", 
                      "Moths", "Odonata", "Orthoptera",
                      "Shieldbugs", "Soldierflies", "Spiders"),
    "dataLink" = c(
      "../Data/Species_data/Raw_data/Butterflies/ccunningham_butterflies.csv", # Butterflies
      "../Data/Species_data/Raw_data/Caddisflies/Caddisflies_2022_Connected Treescapes.rds", # Caddisflies
      "../Data/Species_data/Raw_data/Carabids/Ground_Beetles_2022.rds", # Carabids
      "../Data/Species_data/Raw_data/Centipedes/Centipedes_2022.rds", # Centipedes
      "../Data/Species_data/Raw_data/Ephemeroptera/Ephemeroptera_2022.rds", # Ephemeroptera
      "../Data/Species_data/Raw_data/Gelechiidae/Gelechiidae_2022.rds", # Gelechiidae
      "../Data/Species_data/Raw_data/Hoverflies/Hoverflies_2022.rds", # Hoverflies
      "../Data/Species_data/Raw_data/Ladybirds/Ladybirds_2022.rds", # Ladybirds
      "../Data/Species_data/Raw_data/Molluscs/Molluscs_2022_Connected Treescapes.rds", # Molluscs
      "../Data/Species_data/Raw_data/Moths/ccunningham_moths.csv", # Moths
      "../Data/Species_data/Raw_data/Odonata/Odonata_2022.rds", # Odonata
      "../Data/Species_data/Raw_data/Orthoptera/Orthoptera_2022.rds", # Orthoptera
      "../Data/Species_data/Raw_data/Shieldbugs/Shieldbugs_2022.rds", # Shieldbugs
      "../Data/Species_data/Raw_data/Soldierflies/Soldierflies_2022.rds", # Soldierflies
      "../Data/Species_data/Raw_data/Spiders/Spiders_2022.rds")) # Spiders
              
# CREATE ADDITIONAL COLUMNS AND LOOP THROUGH TAXA ------------------

# Add columns to populate with number of records
taxaFiles$nSpecies <- NA
taxaFiles$nTotalRecords <- NA
taxaFiles$nFilteredRecordsAll <- NA
taxaFiles$nFilteredRecordsPeriod1 <- NA
taxaFiles$nFilteredRecordsPeriod2 <- NA
taxaFiles$nVisitsAll <- NA
taxaFiles$nVisitsPeriod1 <- NA
taxaFiles$nVisitsPeriod2 <- NA

# Loop through each taxaFiles dataframe row
for (i in 1: NROW(taxaFiles)) {
 
  # READ IN DATA ----------------------------------------------------
  
  # Extract file extension
  iExt <- taxaFiles[i, "dataLink"] %>%
    str_sub(., start= -4)

  # Read in data, depending on extension
  if( iExt == ".csv" ) {
    iFile <- taxaFiles[i, "dataLink"] %>%
      read.csv
  }  else if( iExt == ".rds" ) {
    iFile <- taxaFiles[i, "dataLink"] %>%
      readRDS
  } 

  # If taxa groups is molluscs...
  if(taxaFiles[i, "taxaGroups"] == "Molluscs") {
    
    # Carry out additional filtering step
    iFile  <- iFile %>%
      filter(taxon_group_one != "marine molluscs"  | is.na(taxon_group_one)) %>%
      filter(taxon_group_two != "marine molluscs"  | is.na(taxon_group_two))
  }
  
  # COUNT TOTAL ROWS --------------------------------------------------
  
  # Assign row number
  taxaFiles[i, "nTotalRecords"] <- NROW(iFile)
  
  # FILTER ROWS --------------------------------------------------------
  
  # Assign raw data
  taxaData <- iFile
  
  ### STANDARDISE COLUMN NAMES BETWEEN TAXA GROUPS
  
  # Process species name (DEPENDS ON TAXA GROUP)
  if(taxaFiles[i, "taxaGroups"] %in% dataBRC) {
    
    # Remove subgenus "(subgenus)" from a taxon name
    taxaData$recommended_name <- trimws(gsub("\\([[:alpha:]]*\\)[[:space:]]", 
                                             "", 
                                             taxaData$recommended_name))
    # Standardise species name column
    taxaData <- rename(taxaData, 
                       taxon = recommended_name)
  }
  
  ### FILTER RAW DATA
  # N.B. Each row is a record and there can be multiple species records per visit.
  
  # Filter by monad
  
  # Remove records with no monad data
  taxaData <- taxaData[ !is.na(taxaData$monad), ]
  
  # Remove 5x5km records that may have crept into some data requests
  taxaData <- taxaData[grep(pattern = '\\NE$|\\SE$|\\SW$|\\SE$',
                            x = taxaData$monad,
                            invert = TRUE) , ]
  
  # Filter by date
  
  # Format date (DEPENDS ON TAXA GROUP)
  if(taxaFiles[i, "taxaGroups"] %in% dataBC) {
    
    # Remove multi-day records as not usable in this analysis
    taxaData <- taxaData[grep(pattern = '-',
                              x = taxaData$date,
                              invert = TRUE) ,]
    # Reformat date
    taxaData$date <- as.POSIXct(taxaData$date, format = "%d/%m/%Y")
  }
  
  if(taxaFiles[i, "taxaGroups"] %in% dataBRC) {
    
    # Remove multi-day records as not usable in this analysis
    taxaData <- taxaData[taxaData$startdate == taxaData$enddate, ]
    
    # Create new standardised "date" column or, if column already exists, overwrite
    taxaData <- mutate(taxaData, date = startdate)
    
    # Reformat date
    taxaData$date <- as.POSIXct(taxaData$date, format = "%Y-%m-%d")
  }
  
  ## Add year column
  taxaData$year <- taxaData$date %>%
    format(., format="%Y") %>% 
    as.numeric(.)
  
  # Aggregate annual data to multi-year groups
  taxaData$iYear <- sapply(taxaData$year, function(x) {
    
    assignYearGroup(x,
                    yearGroupRange1 = range1,
                    yearGroupRange2 = range2) })
  
  # Drop years outside year groups
  taxaData <- taxaData[!is.na(taxaData$iYear),]
  
  # Select useful columns
  taxaData <- dplyr::select(taxaData, taxon, monad, date, year, iYear)
  
  # COUNT FILTERED RECORDS AND SPECIES ------------------------------
  
  # Assign records count
  taxaFiles[i, "nFilteredRecordsAll"] <- NROW(taxaData)
  taxaFiles[i, "nFilteredRecordsPeriod1"] <- NROW(taxaData[taxaData$iYear == 1,])
  taxaFiles[i, "nFilteredRecordsPeriod2"] <- NROW(taxaData[taxaData$iYear == 2,])

  # Assign species/aggregates number
  taxaFiles[i, "nSpecies"] <- unique(taxaData$taxon) %>%
    length
  
  # PROCESS TO VISITS DATAFRAME ---------------------------------
  
  ### GET SPECIES NAME
  
  # Set species
  # N.B. Here we are just counting visits and not interested in which species,
  # so species can be any. Here we arbitrarily set to 1
  iSpecies <- unique(taxaData$taxon)[1]
  
  ### ADD DATA COLUMNS
  
  # Add visit column (concatenate gridref plus date)
  taxaData$visit <- paste(taxaData$monad, taxaData$date, sep = "_")
  
  # Add iSpecies presence column
  taxaData$presence <- if_else(taxaData$taxon == iSpecies, 1, 0)
  
  # Drop taxon column
  taxaData <- dplyr::select(taxaData, -taxon)
  
  # Add count of number of records at each visit
  # (Number of records reported on a visit, as an indicator of effort -
  # most records are of list length 1 and probably not a comprehensive survey)
  taxaData <- add_count(taxaData, visit, name = "numRecords")
  
  ### REDUCE ROW COUNT TO UNIQUE VISITS
  
  # Filter all records to a single value for each visit
  visitData <- taxaData %>%
    # For each visit ...
    group_by(visit) %>%
    # Extract max presence value (i.e. if any records are iSpecies then 1, else 0)
    slice(which.max(presence)) %>%
    # Finally, ungroup
    ungroup()

  # COUNT VISITS --------------------------------------------------
  
  # Assign row numbers
  taxaFiles[i, "nVisitsAll"] <- NROW(visitData)
  taxaFiles[i, "nVisitsPeriod1"] <- NROW(visitData[visitData$iYear == 1,])
  taxaFiles[i, "nVisitsPeriod2"] <- NROW(visitData[visitData$iYear == 2,])
  
  # CLEAR MEMORY --------------------------------------------------
  
  # Remove files and clear memory
  rm(iFile, taxaData, visitData)
  gc()
  
}

# CALCULATE SUMMARY STATISTICS ------------------------------------

# Save summary data
save(taxaFiles, file = "../Data/Species_data/taxaRecordSummary.RData")

# Sums across taxa
colSums(taxaFiles[,-c(1:2)])


