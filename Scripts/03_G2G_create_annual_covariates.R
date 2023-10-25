# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Calculate an annual soil moisture metric from G2G data
#
# Script Description: Using the downloaded G2G data, process the 
# soil moisture data by taking annual means of soil moisture
# 
# Available here: 
# https://catalogue.ceh.ac.uk/documents/c9a85f7c-45e2-4201-af82-4c833b3f2c5f

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to C: (permission issues on R:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)
library(ncdf4)

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/bng.prj")

bng <- sf::st_crs("../Data/Spatial_data/bng.prj")$wkt

### LIST DIRECTORIES ------------------------------

###  G2G data directory
G2GDir <- "../Data/Spatial_data/G2G/Observed/data"

# Set GB file name
fileGB <- paste0(G2GDir, "/g2g_gb_mmsoil_obs_1980_2011.nc")

# Set GB lakes file name
fileLakesGB <- paste0(G2GDir, "/ukscape_g2g_gb_soilmoisture_lakegrid.nc")
 
# Set NI file name
fileNI <- paste0(G2GDir, "/g2g_ni_mmsoil_obs_1980_2011.nc")

# Set NI lakes file name
fileLakesNI <- paste0(G2GDir, "/ukscape_g2g_ni_soilmoisture_lakegrid.nc")

### Output directories

# Output directory for GDD5
destDir <- "../Data/Spatial_data/Annual_climate_covar/soilM"

# Create output directories if they do not exist
if(!file.exists(destDir)) {dir.create(destDir) }

### READ G2G (SOIL MOISTURE) DATA ---------------------------- 
# N.B. Time is stored as time bounds (i.e. start and end of month). 
# However, the starting time bound is the day before the first of the 
# month, and this is what is automatically  extracted by rast(). This is
# unhelpful for filtering years, months etc. Hence the ending time bound
# has to be manually extracted from the .nc file and added.

### GB DATA

### Create rasters

# Create G2G spatRaster from .nc file
spatRastGB <- rast(fileGB)

# Create majority lake cover spatRaster from .nc file
spatRastLakesGB <- rast(fileLakesGB)

# Assign CRS
crs(spatRastGB) <- crs(spatRastLakesGB) <- bng

### Change time to ending time bound

# Read .nc file
ncGB <- nc_open(fileGB)

# Extract ending time bound (this is more helpful as within month)
timeGB <- ncvar_get(ncGB,"time_bnds")[2,] %>%
  as.Date(., origin = '1900-01-01') # Convert to year-month-day format

# Manually change layer names to ending time bound
names(spatRastGB) <- as.character(timeGB)

### NI DATA

### Create rasters

# Create G2G spatRaster from .nc file
spatRastNI <- rast(fileNI)

# Create majority lake cover spatRaster from .nc file
spatRastLakesNI <- rast(fileLakesNI)

# Assign CRS
crs(spatRastNI) <- crs(spatRastLakesNI) <- bng

### Change time to ending time bound

# Read .nc file
ncNI <- nc_open(fileNI)

# Extract ending time bound (this is more helpful as within month)
timeNI <- ncvar_get(ncNI,"time_bnds")[2,] %>%
  as.Date(., origin = '1900-01-01') # Convert to year-month-day format

# Manually change layer names to ending time bound
names(spatRastNI) <- as.character(timeNI)

### PROCESS G2G (SOIL MOISTURE) DATA ---------------------------- 

# REMOVE INCOMPLETE YEARS

# Check whether all layer times are identical for GB and NI,
# i.e. can be correctly compared and merged...
if(identical(names(spatRastGB),names(spatRastNI))) {
  
  # ..if they are identical: find which years are incomplete
  yearsExclude <- names(spatRastGB) %>% 
    substr(.,1,4) %>% # Extract just year from layer names
    table(.)%>% # Convert to table format
    '<' (12) %>% # Find years which have less than 12 layers
    which(.) %>%
    names(.)
  
  } else {
    
    # ...if not:
    print("Different times for GB and NI")
  }

# Set string to use to remove incomplete years
yearsExcludeStr <- paste0( yearsExclude, "-")

# Find spatRast layers to include, i.e. inverse of yearsExclude
includeLayers <- grep( paste(yearsExcludeStr, collapse = "|"),
                       names(spatRastGB),
                       invert = TRUE)

# Subset GB and NI SpatRasters
spatRastGB <- spatRastGB[[includeLayers]]
spatRastNI <- spatRastNI[[includeLayers]]

### MERGE

# Merge G2G spatRasters together
spatRastUK <- mosaic(spatRastGB, spatRastNI)

# Merge lake spatRasters
spatRastLakesUK <- mosaic(spatRastLakesGB, spatRastLakesNI)

### PROCESS LAKES SPATRASTER

# Change extent of lakes spatRaster to G2G
spatRastLakesUK <- crop(spatRastLakesUK, spatRastUK)

# Change from 2s and 1s, to just 1s and NAs
spatRastLakesUK <- ifel(spatRastLakesUK == 2,
                        yes = 1,
                        no = NA)

# Save lakes spatRaster
writeRaster(spatRastLakesUK,
            "../Data/Spatial_data/G2G/Lakes.tif",
            overwrite = TRUE)

# CREATE ANNUAL RASTERS --------------------------------------- 

# Find all complete years
years <- names(spatRastUK) %>%
  substr(.,1,4) %>%
  unique(.)

# Loop through each year and produce annual variable separately
# N.B. Here we simply take the annual mean, including every month
for (i in years) {

  # Subset to months (layers) within year
  spatRastUKi <- spatRastUK[[ grep( i,
                                    substr( names(spatRastUK), 1, 4)) ]]
  
  # Calculate annual mean
  soilMean <- mean(spatRastUKi)
  
  # Write raster
  writeRaster(soilMean,
              paste0(destDir, "/", "soilM_", i, ".tif"),
              overwrite = TRUE)
}
