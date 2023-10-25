# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Species record comparison by country
#
# Script Description: Used to compare recording effort in NI
# to GB as scoping exercise
#

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)

# DATA FILES ------------------------------------------

# GB data
recordsGB <- read.table("../Data/Species_data/taxaSummaryGB.txt", header = TRUE, sep = '\t')  

# NI data 
recordsNI <- read.table("../Data/Species_data/taxaSummaryNI.txt", header = TRUE, sep = '\t')

# JOIN DATA FILES TOGETHER ------------------------------

# Add Country column to recordsNI
recordsNI$COUNTRY <- "Northern Ireland"

# Remove NI from GB dataset
recordsGB <- recordsGB[recordsGB$COUNTRY %in% c("England", "Scotland", "Wales"),] 

# Join together
recordsUK <- full_join(recordsGB, recordsNI, by = c("Taxa", "nuRecs", "nu1kmGrids", "COUNTRY"))

# CALCULATE RECORDS AND COVERAGE VALUES -----------------

# Country approx sizes
EnglandKm2 <- 130000; ScotlandKm2 <-78000; WalesKm2 <- 21000; NIKm2 <- 14000; 

### PROPORTION COVERAGE

# Create coverage column
recordsUK$percent1kmGrids <- NA

# Fill Column
recordsUK$percent1kmGrids <- sapply(1:NROW(recordsUK), FUN = function(x) {
  
  if(recordsUK[x,"COUNTRY"] == "England") { 100* recordsUK[x,"nu1kmGrids"] / EnglandKm2  }
  
  else if (recordsUK[x,"COUNTRY"] == "Scotland") { 100 * recordsUK[x,"nu1kmGrids"] / ScotlandKm2 }
  
  else if (recordsUK[x,"COUNTRY"] == "Wales") { 100 * recordsUK[x,"nu1kmGrids"] / WalesKm2 } 
  
  else if (recordsUK[x,"COUNTRY"] == "Northern Ireland") { 100* recordsUK[x,"nu1kmGrids"] / NIKm2 } 
})


### RECORDS PER AREA

# Create coverage column
recordsUK$recordsToArea <- NA

# Fill Column
recordsUK$recordsToArea <- sapply(1:NROW(recordsUK), FUN = function(x) {
  
  if(recordsUK[x,"COUNTRY"] == "England") { recordsUK[x,"nuRecs"] / EnglandKm2  }
  
  else if (recordsUK[x,"COUNTRY"] == "Scotland") { recordsUK[x,"nuRecs"] / ScotlandKm2 }
  
  else if (recordsUK[x,"COUNTRY"] == "Wales") { recordsUK[x,"nuRecs"] / WalesKm2 } 
  
  else if (recordsUK[x,"COUNTRY"] == "Northern Ireland") { recordsUK[x,"nuRecs"] / NIKm2 } 
})

### VISUALISE ---------------------------------------

recordsUKusable <- recordsUK[!recordsUK$Taxa %in% 
                               c("AquaticBugs", "Fish", "LeafSeedBeetles", "PlantBugs", "Weevils", "SoldierBeetles", "RoveBeetles",
                                 "HypogeanCrustacea", "FungusGnats", "E&D", "Plecoptera", "Neuropterida" ),]

# Coverage
ggplot() +
  geom_bar(data = recordsUKusable, aes(x = Taxa, y = percent1kmGrids, fill = COUNTRY),
           position = "dodge", stat = "identity")  

# Records
ggplot() +
  geom_bar(data = recordsUKusable, aes(x = Taxa, y = sqrt(recordsToArea), fill = COUNTRY),
           position = "dodge", stat = "identity")  

