# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: [CLUSTER SCRIPT]
# Calculate the three temperature metrics from CHESS met data
#
# Script Description: Using the downloaded CHESS met data, process the 
# daily surface air temperature .nc files into usable raster stacks and calculate
# annual metrics. These are minimum winter temperature, growing degree days and seasonality.
# N.B. temperature units in Kelvin (K)

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to C: (permission issues on R:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)

### SET PARAMETERS AND LIST DIRECTORIES ------------------------------

# Each iteration is one year, specified array in job script
args <- commandArgs(trailingOnly = TRUE)
iYear <- as.numeric(args[1])

### Species identify starts here:

### CHESS-met temperature (tas) data

# Specify CHESS folder
chessDir <- "Treescapes/Climate_processing/CHESS"

# Specify all folders where monthly data stored 
# (2018 and 2019 are seperate as had to be downloaded manually)
tasFolder1970_2017 <- paste0( chessDir, "/CHESS_met/Current/chess_tas_1970_2017")
tasFolder2018 <- paste0( chessDir, "/CHESS_met/Current/chess_tas_2018")
tasFolder2019 <- paste0( chessDir, "/CHESS_met/Current/chess_tas_2019")

### Output directories

# Output directory for GDD5
destDirGDD5 <- paste0( chessDir, "/Current_climate_variables/GDD5")

# Output directory for MTCO
destDirMTCO <- paste0( chessDir, "/Current_climate_variables/MTCO")

# Output directory for tasCV
destDirTasCV <- paste0( chessDir, "/Current_climate_variables/tasCV")

# Create output directories if they do not exist
lapply ( c(destDirGDD5, destDirMTCO, destDirTasCV), function(x) { 
          
  if(!file.exists(x)) {
    dir.create(x)
    }
})

### LOAD CHESS-MET TEMPERATURE (TAS) DATA ------------------------------

# Note each .nc file contains summarised daily values for each month, 
# i.e. a stack of around 30 daily values for every 1x1km cell in UK

# List all .nc files within folders
tasFiles1970_2017 <- list.files(tasFolder1970_2017, full.names = TRUE)
tasFiles2018 <- list.files(tasFolder2018, full.names = TRUE)
tasFiles2019 <- list.files(tasFolder2019, full.names = TRUE)

# Bind all file locations together
tasFilesAll <- c(tasFiles1970_2017, tasFiles2018, tasFiles2019)

# Create year associated with each .nc file 
tasFilesYear <- sub('.*\\/', '', tasFilesAll) %>% # Drop text before last "/"
  sub("chess-met_tas_gb_1km_daily_", "", .) %>% # Drop text before "year"
  substr(., 1,4) %>% # Extract the 4 'year' characters
  as.integer(.) # Convert to number
  
# Create list of all years
years <- unique(tasFilesYear)

### Read in .nc files

# Subset files to only year 'i'
tasFiles <- tasFilesAll[tasFilesYear == years[iYear]]
  
# Rasterise all files into same raster brick
tasRast <- rast(tasFiles)

### PROCESS TEMPERATURE ENVIRONMENTAL VARIABLES ---------------------------

### GROWING DEGREE DAYS (GDD5)

# Measure of heat accumulation, predictor of development rates
# Calculated as the sum of daily temperatures greater than 5 Celsius

# Eqn:
# (1) daily GDD = tas - Tbase [Tbase set to 5 celsius]
# Note negative values ignored
# (2) annual GDD = sum(daily GDD)

# Convert 5 degrees Celsius to Kelvin (CHESS unit)
fiveCinK <- 5 + 273.15 

# Subtract 5 degrees (in Kelvin)
allyear.m5 <- tasRast - fiveCinK

# Set days below 5 degrees to 0 so they will be ignored
allyear.m5[allyear.m5 < 0] <- 0

# Sum
GDD5 <- sum(allyear.m5) # sum temp above 5oC

# Write raster
writeRaster(GDD5,
            paste0(destDirGDD5, "/", "GDD5_", years[iYear], ".tif"),
            overwrite = TRUE)

### MEAN DAILY TEMPERATURE OF THE COLDEST MONTH (MTCO)

# First, loop through monthly files of daily data to find the monthly temperature means
tasMonthlyMean <- lapply(tasFiles, function(x) {
  
  # Rasterize
  tasMonth <- rast(x)
  
  # Find mean, and return value
  mean(tasMonth) %>%
    return(.)

}) %>%
  rast(.) # Combine output list into single spatRaster

# Find the coldest monthly mean temperature for each cell
MTCO <- min(tasMonthlyMean)

# Write raster
writeRaster(MTCO,
            paste0(destDirMTCO, "/", "MTCO_", years[iYear], ".tif"),
            overwrite = TRUE)

### SEASONALITY

# Get cell-wise standard deviation
tasVar <- app(tasRast, sd) 

# Get cell-wise mean
tasMean <- mean(tasRast)

# Get annual seasonality coefficient (sd / mean)
tasCV <- tasVar / tasMean  

# Write raster
writeRaster(tasCV,
            paste0(destDirTasCV, "/", "tasCV_", years[iYear], ".tif"),
            overwrite = TRUE)
