# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: 
#
# Script Description: Need to:
#
# Read in Land cover data 25m for 1990 and 2020


# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Install.packages
# install.packages("janitor")

# Load packages
library(terra)
library(tidyverse)

# DATA FILES ------------------------------------------

# Land cover map folder
dirLCM <- "../Data/Spatial_data/LCM/"

# 1990 25m Land Cover Map
LCM1990File <- "1990_25m/gb1990lcm25m.tif"

# 2020 25m Land Cover Map
LCM2020File <-"2020_25m/gb2020lcm25m.tif"

# Study area boundary polygon

testAreaFile <- "../Data/Spatial_data/NFC boundary/Forest boundary.shp"

# READ DATA ------------------------------------------

# Read LCM rasters
LCM1990all <- rast(paste0(dirLCM, LCM1990File))
LCM2020all <- rast(paste0(dirLCM, LCM2020File))

# Read study area
studyArea <- vect(testAreaFile)

# Limit to study area for testing
LCM1990 <- terra::crop(LCM1990all, studyArea)
LCM2020 <- terra::crop(LCM2020all, studyArea)


# TREESCAPE COVER ------------------------------------

# Get broadleaf cover
LCM1990BF <- LCM1990[[1]] == 1
LCM2020BF <- LCM2020[[1]] == 1

# Get coniferous cover
LCM1990CF <- LCM1990[[1]] == 2
LCM2020CF <- LCM2020[[1]] == 2

# Aggregate to 1km

test <- aggregate(LCM1990BF, 40, sum)
plot(test)


LCM1990BF %>% plot(.)










