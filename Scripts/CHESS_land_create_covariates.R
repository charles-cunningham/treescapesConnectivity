# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Calculate the soil moisture metric from CHESS land data
#
# Script Description: Using the downloaded CHESS land data, process the 
# soil moisture data.

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to C: (permission issues on R:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)
library(ncdf4)

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.prettywkt?download",
              destfile = "../Data/Spatial_data/bng.prj")

bng <- "../Data/Spatial_data/bng.prj"

### SET PARAMETERS AND LIST DIRECTORIES ------------------------------

# Each iteration is one year
iYear <- 1

###  Soil moisture (HydEn) data
landDir <- "../Data/Spatial_data/CHESS/CHESS_land/Current"

### Output directories

# Output directory for GDD5
destDirSoil <- "../Data/Spatial_data/CHESS/Current_climate_variables/SoilM"

# Create output directories if they do not exist
if(!file.exists(destDirSoil)) {dir.create(destDirSoil) }

### FIND CHESS-LAND MOISTURE (HYDEN) DATA ---------------------------- 

# Specify all files where monthly data stored 
# (2016 to 2019 are separate as had to be downloaded manually)
landFiles <- list.files(landDir,
                        pattern = "\\.nc$",
                        recursive = TRUE)

# Create year associated with each .nc file 
soilFilesYear <- sub("\\.nc$", "", landFiles)  %>% # Drop .nc
  substr(., nchar(.) - 3 ,nchar(.)) %>% # Extract last 4 'year' characters
  as.integer(.) # Convert to number

# Create list of all years
years <- unique(soilFilesYear)

### LOAD CHESS-LAND MOISTURE (HYDEN) DATA ---------------------------

# Begin loop through years at this point
for (i in years) {

# Subset files to only year 'i'
soilNC <- paste0(landDir, "/", landFiles[soilFilesYear == i])

# Read in .nc file
soilData <- nc_open(soilNC)

# Extract the variable of interest ( using Gridbox soil moisture availability factor (beta) [fsmc_gb] for now)
soilVar <- ncvar_get(soilData,"fsmc_gb")

# Extract x and y
x <- ncvar_get(soilData, "eastings" )
y <- ncvar_get(soilData, "northings")

# Create a dataframe
soilDF <- data.frame(x,y,soilVar)

# Rasterise
soilR <- rast(soilDF, type="xyz", crs=bng)

# Continue loop to next section

### PROCESS SOIL MOISTURE ENVIRONMENTAL VARIABLES ---------------------------

# Calculate annual mean
soilMean <- mean(soilR)

# Write raster
writeRaster(soilMean,
            paste0(destDirSoil, "/", "soilM_", i, ".tif"),
            overwrite = TRUE)

}
