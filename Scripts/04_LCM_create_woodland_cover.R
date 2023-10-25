# HEADER -------------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Create woodland cover maps
#
# Script Description: Create woodland cover map for use in modelling -
# UK 1x1km resolution using CEH Land Cover Map (LCM) products.
# Approach is as follows, within loop:
#
# - Aggregate 1990 25m LCM to 1x1km for broadleaf, coniferous and both ('woodland')
# - Aggregate 2015 25m LCM to 1x1km for broadleaf, coniferous and both ('woodland')
# - Do same for both 2020 10m and 25m LCM.
# - Use 2020 10m and 25m to understand relationship between the spatial resolutions

# N.B. We only have LCM data for 1990 and 2015 that are comparable (unfortunately 2000
# and 2007 LCM are not comparable due to different methodologies).
# N.B.B.  BF = broadleaf woodland, CF = coniferous woodland, W = all woodland

# LOAD LIBRARIES & INSTALL PACKAGES ----------------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(tidyverse)
library(ggplot2)

# DATA FILES ---------------------------------------------

# Land cover map folder
dirLCM <- "../Data/Spatial_data/LCM/"

### 1990 Land Cover Map

# 1990 25m GB (no 10m equivalent)
GB1990File <- "1990/GB_25m/gb1990lcm25m.tif"
# 1990 25m NI
NI1990File <- "1990/NI_25m/LCM1990_NI_25m.tif"

### 2015 Land Cover Map

# 2015 25m GB (no 10m equivalent)
GB2015File <- "2015/GB_25m/lcm2015gb25m.tif"
# 2015 25m NI
NI2015File <- "2015/NI_25m/lcm2015ni25mfinal.tif"

### 2020 Land Cover Map

# 2020 25m GB
GB25m2020File <-"2020/GB_25m/gb2020lcm25m.tif"
# 2020 25m NI
NI25m2020File <-"2020/NI_25m/nilcm25m2020.tif"
# 2020 10m GB
GB10m2020File <-"2020/GB_10m/gb2020lcm10m.tif"
# 2020 10m NI
NI10m2020File <-"2020/NI_10m/ni2020lcm10m.tif"

### CODE --------------------------------------------------

# WOODLAND LAND COVER MAPS FOR 1990 2015, and 2020 --------
  
# Create table of files, year, and resolution to populate loop
LCMdf <- data.frame("year"   = c("1990",        "2015",       "2020",          "2020"),
                    "res"    = c("25",          "25",         "25",            "10"),
                    "GBfile" = c("GB1990File" , "GB2015File", "GB25m2020File", "GB10m2020File"),
                    "NIfile" = c("NI1990File" , "NI2015File", "NI25m2020File", "NI10m2020File"))

# Start loop here
# (complete loop takes approximately 4 hours to run if all files are created)
for (i in 1:NROW(LCMdf)) {

  # Assign year and resolution for row i
  year <- LCMdf[i,"year"]
  res <- LCMdf[i,"res"]

  # First, check if .tifs created for each year in this section already exist..
  if(  !file.exists( paste0(dirLCM, year, "/",
                           "UK", year, "all_1km.tif")) &
       !file.exists( paste0(dirLCM, year, "/",
                            "UK", res, "m", year, "all_1km.tif"))) {
      
    # Read in GB spatRast
    GB <-  get(LCMdf[i, "GBfile"]) %>% # Get full file name
      paste0(dirLCM, .) %>% # Get full file path
      rast(.) %>% # Create spatRast
      .[[1]] # Band 1 is LCM type
     
    # Read in NI spatRast
    NI <-  get(LCMdf[i, "NIfile"]) %>% # Get full file name
      paste0(dirLCM, .) %>% # Get full file path
      rast(.) %>% # Create spatRast
      .[[1]] # Band 1 is LCM type    
               
    # Change the CRS of NI1990 from OSNI to OSGB using GB1990 as template
    NI <- project(NI, GB)

    # Re-classify 0 to NA for GB and NI before merging
    GB <- classify(GB, cbind( 0, NA ))
    NI <- classify(NI, cbind( 0, NA ))

    # Merge GB and NI LCM into single UK LCM
    UK <- mosaic(GB, NI)

    # Get broadleaf cover
    UKbf <- classify(UK,
                     cbind( 1, 1 ), # Reassign broadleaf [value = 1] to '1'
                     others = 0) # and other classes to 0

    # Get coniferous cover
    UKcf <- classify(UK,
                     cbind( 2, 1 ), # Reassign coniferous [value = 2] to '1'
                     others = 0) # and other classes to 0

    # Get total woodland cover
    UKw <- classify(UK,
                    cbind( c(1,2), 1 ), # Reassign deciduous and coniferous to '1'
                    others = 0) # and other classes to 0

    # Stack to single raster
    UKall <- c(UKbf, UKcf, UKw)

    # Set names
    names(UKall) <- c("BF", "CF", "W")
  
    # Write spatRast as a single multi-band raster
    if (year == "2020") { # Include resolution in file name for 2020
      
      writeRaster(UKall,
                  paste0( dirLCM, year, "/",
                          "UK", res, "m", year, "all.tif"),
                  overwrite = TRUE)
      } else {
        
        writeRaster(UKall,
                    paste0( dirLCM, year, "/",
                            "UK", year, "all.tif"),
                    overwrite = TRUE)
      }
    
    # Calculate the aggregation factor from current resolution to 1x1km cells
    aggFact <- 1000 / as.integer( res )
    
    # Aggregate all to 1x1km scale through proportion cell coverage
    UKall_1km <-  aggregate(UKall, aggFact, sum, na.rm = TRUE) / (aggFact^2)
    
    # Write 1x1km aggregated spatRast as a single multi-band raster
    if (year == "2020") { # Include resolution in file name for 2020
      
      writeRaster(UKall_1km,
                  paste0( dirLCM, year, "/",
                          "UK", res, "m", year, "all_1km.tif"),
                  overwrite = TRUE) 
    } else {
      
      writeRaster(UKall_1km,
                  paste0( dirLCM, year, "/",
                          "UK", year, "all_1km.tif"),
                  overwrite = TRUE) 
    }
    
    # Remove spatRasts, and free up memory
    rm(GB, NI, UK,
       UKbf, UKcf, UKw,
       UKall, aggFact, UKall_1km)
    
    gc()
  
}}

# INVESTIGATE 2020 10M AND 25M RASTER

### Load in 1x1km rasters

# 2020 25m
UK25m2020all_1km <- rast( paste0( dirLCM, "2020/UK25m2020all_1km.tif"))
## 2020 10m
UK10m2020all_1km <- rast( paste0( dirLCM, "2020/UK10m2020all_1km.tif"))

# Change names of 10m spatRast to differentiate between the two
names(UK10m2020all_1km) <- c("BF_10m", "CF_10m", "W_10m")

# Create data.frame of 25m and 10m coverage
compare10m25m <- data.frame(UK25m2020all_1km[], UK10m2020all_1km[]) %>% na.omit(.)

# Compare broadleaf
ggplot( data = compare10m25m, aes( x = compare10m25m$BF, y = compare10m25m$BF_10m)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) + 
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme_classic() +
  theme()

# Compare coniferous
ggplot( data = compare10m25m, aes( x = compare10m25m$CF, y = compare10m25m$CF_10m)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) + 
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme_classic() +
  theme()

# Compare all woodland
ggplot( data = compare10m25m, aes( x = compare10m25m$W, y = compare10m25m$W_10m)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE) + 
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  theme_classic() +
  theme()
