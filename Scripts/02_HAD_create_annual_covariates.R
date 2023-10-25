# HEADER ---------------------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: [CLUSTER SCRIPT]
# Calculate climate metrics from HAD data
#
# Script Description: Using the downloaded HAD data, process the monthly rainfall, 
# daily minimum and maximum surface air temperature .nc files into usable spatRasters.
# Calculate annual metrics from these: total annual precipitation (RAIN), minimum winter 
# temperature (MTCO), growing degree days (GDD5) and seasonality (tasCV). 
# N.B. temperature units in Kelvin (K) for seasonality, otherwise Celsius.

### LOAD LIBRARIES & INSTALL PACKAGES ---------------------------------

# Change  library to R: (not enough space on C:)
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(tidyverse)

# Set fraction of RAM that may be used by the program
terraOptions(memfrac = 0.9)

### SET PARAMETERS AND LIST DIRECTORIES ------------------------------

# Each iteration is one year, specified array in job script
args <- commandArgs(trailingOnly = TRUE)
iYear <- as.numeric(args[1])

### Specify input directories

# Data directory
dataDir <- "Treescapes/Climate_processing"

# Individual variable sub-directories
tasminDir <-  paste0( dataDir, "/HAD/tasmin") # Daily minimum temperature
tasmaxDir <-  paste0( dataDir, "/HAD/tasmax") # Daily maximum temperature
rainfallDir <- paste0( dataDir, "/HAD/rainfall") # Monthly rainfall

### Specify output directories

# Output directory for GDD5
destDirGDD5 <- paste0( dataDir, "/Annual_climate_covar/GDD5")

# Output directory for WMIN
destDirWMIN <- paste0( dataDir, "/Annual_climate_covar/WMIN")

# Output directory for tasCV
destDirTasCV <- paste0( dataDir, "/Annual_climate_covar/tasCV")

# Output directory for RAIN
destDirRAIN <- paste0( dataDir, "/Annual_climate_covar/RAIN")

# Create output directories if they do not exist
lapply (c(destDirGDD5, destDirWMIN, destDirTasCV, destDirRAIN),
        function(x) {
          
          if (!file.exists(x)) {
            dir.create(x, recursive = TRUE) }
      })

### Download BNG WKT string
# N.B. Individual filename needed for each task to prevent
# different tasks trying to read/write at same time 

# Specify temporary file to download 'bng' CRS wkt to
tempFile <- paste0("Treescapes/Climate_processing/",iYear, "bng.prj")

# Download file
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = tempFile)

# Assign wkt string
bng <- sf::st_crs(tempFile)$wkt

# Remove temporary file
unlink(tempFile)

### LOAD ENVIRONMENTAL VARIABLES ----------------------------------------------

### LOAD MINIMUM AND MAXIMUM TEMPERATURE DATA 
# Note each .nc file contains summarised daily values for each month, 
# i.e. a stack of around 30 daily values for every 1x1km cell in UK

### LIST FILES

# List all .nc files within folders
tasminFiles <- list.files(tasminDir, full.names = TRUE)
tasmaxFiles <- list.files(tasmaxDir, full.names = TRUE)

# Create year associated with tasmin .nc files
tasminFilesYear <- str_extract(tasminFiles, # Extract 'year' characters
                               "(?<=tasmin_hadukgrid_uk_1km_day_).{4}") %>% 
  as.integer(.) # Convert to number

# Create year associated with tasmax .nc files
tasmaxFilesYear <- str_extract(tasmaxFiles, # Extract 'year' characters
                               "(?<=tasmax_hadukgrid_uk_1km_day_).{4}") %>% 
  as.integer(.) # Convert to number

# Check if tasmax and tasmin file years are identical, and...
if(identical(tasminFilesYear, tasmaxFilesYear)) {
  
  # ...if they are, then create list of all years
  years <- unique(tasminFilesYear)
  
} else {
  
  # ...if not:
  print("Different years for tasmin and tasmax")
}

# Find number of years, this will be the ending array value for the job script
length(years)
# N.B. Starting array value is '2' as t-1 year is needed for WMIN calculation,
# i.e. array = 2:43. Hence why we need 1979 to start from 1980

### READ IN IYEAR .NC FILES

# Subset files to only year 'i', and rasterise
tasminR <- tasminFiles[tasminFilesYear == years[iYear]] %>%
  rast(.)
tasmaxR <- tasmaxFiles[tasmaxFilesYear == years[iYear]] %>%
  rast(.)  

# Assign CRS
crs(tasminR) <- crs(tasmaxR) <- bng

# Change names to date
names(tasminR) <- sub(" .*", "", time(tasminR))
names(tasmaxR) <- sub(" .*", "", time(tasmaxR))

### READ IN DECEMBER OF (IYEAR - 1) .NC FILE
# N.B. We need tasmin December of the previous year for winter cold calculation

# Subset only December of year 'i-1' file, and rasterise
tasminLastDecR <- tasminFiles[tasminFilesYear == years[iYear-1]][12] %>%
  rast(.)

# Assign CRS
crs(tasminLastDecR) <- bng

# Change names to date
names(tasminLastDecR) <- sub(" .*", "", time(tasminLastDecR))

### PROCESS DAILY MEAN TEMPERATURE (TAS)

# Calculate iYear daily mean temperature (only needed for iYear)
tasR <- mean(tasminR, tasmaxR)

### LOAD RAINFALL DATA
# Note each .nc file contains summarised monthly values for each year, 
# i.e. a stack of 12 monthly values for every 1x1km cell in UK

# List all .nc files within folders
rainfallFiles <- list.files(rainfallDir, full.names = TRUE)

# Create year associated with rainfall .nc files
rainfallFilesYear <- str_extract(rainfallFiles, # Extract 'year' characters
                                 "(?<=rainfall_hadukgrid_uk_1km_mon_).{4}") %>% 
  as.integer(.) # Convert to number

# Subset files to only year 'i', and rasterise
rainfallR <- rainfallFiles[rainfallFilesYear == years[iYear]] %>%
  rast(.)

# Assign CRS
crs(rainfallR) <-  bng

# Change names to date
names(rainfallR) <- sub(" .*", "", time(rainfallR))

### PROCESS ENVIRONMENTAL VARIABLES -------------------------------

### TOTAL ANNUAL RAINFALL

# Measure of annual precipitation only. Does not account for topology,
# evaporation, geology as soil moisture does

# Get cell-wise sum
rainfallSum <- sum(rainfallR)

# Write raster
writeRaster(rainfallSum,
            paste0(destDirRAIN, "/", "RAIN_", years[iYear], ".tif"),
            overwrite = TRUE)

### GROWING DEGREE DAYS (GDD5)

# Measure of heat accumulation, predictor of development rates
# Calculated as the sum of daily temperatures greater than 5 Celsius

# Eqn:
# (1) daily GDD = tas - Tbase [Tbase set to 5 celsius]
# Note negative values ignored
# (2) annual GDD = sum(daily GDD)

# Subtract 5 degrees
allyear.m5 <- tasR - 5

# Set days below 5 degrees to 0 so they will be ignored
allyear.m5 <- clamp(allyear.m5, lower = 0)

# Sum temp above 5oC
GDD5 <- sum(allyear.m5)

# Write raster
writeRaster(GDD5,
            paste0(destDirGDD5, "/", "GDD5_", years[iYear], ".tif"),
            overwrite = TRUE)

### MEAN MINIMUM DAILY TEMPERATURE OF THE COLDEST MONTH (WMIN)

# For this measure we use December of the previous year, and 'tasmin' not 'tas',
# as WMIN is a measure of 'winter' cold.

# Hence from year i we take months 1-11, and 12 from previous year
iYearMonthsWMIN <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11")

# Loop through iYear months to find the monthly mean *minimum* temperature
tasminMonthlyMean <- lapply(iYearMonthsWMIN, function(x) {
  
  # For month 'x', extract associated daily tasmin spatRast layers
  tasminMonth <- tasminR[[ grep( paste0("-", x, "-"),
                                 names(tasminR) )  ]]
  
  # Find mean, and return value
  mean(tasminMonth) %>%
    return(.)
  
}) %>%
  rast(.) # Combine output list into single spatRaster

# Find mean of previous December
tasminLastDecMean <- mean(tasminLastDecR)

# Combine year i monthly means with previous December mean
tasminMonthlyMeanAll <- c(tasminMonthlyMean, tasminLastDecMean)

# Find the coldest monthly mean *minimum* temperature for each cell
WMIN <- min(tasminMonthlyMeanAll)

# Write raster
writeRaster(WMIN,
            paste0(destDirWMIN, "/", "WMIN_", years[iYear], ".tif"),
            overwrite = TRUE)

### SEASONALITY

# Convert to Kelvin for seasonality calculation
tasR_K <- tasR + 273.15 

# Get cell-wise standard deviation
tasVar <- app(tasR_K, sd)

# Get cell-wise mean
tasMean <- mean(tasR_K)

# Get annual seasonality coefficient (sd / mean)
tasCV <- tasVar / tasMean  

# Write raster
writeRaster(tasCV,
            paste0(destDirTasCV, "/", "tasCV_", years[iYear], ".tif"),
            overwrite = TRUE)
