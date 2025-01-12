# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Process spatial covariates for inlabru modelling 
#
# Script Description: This script processes the following for downstream 
# use in the analysis:
# - boundaries for masking species records and modelling domain, predictions etc.
# - woodland cover data (previously processed LCM data)
# - connectivity data (output from Omniscape run)
# - annual climate covariate data

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to R: (C: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra) # terra' must be loaded before 'sf' due to 'geodata' dependency bug
library(sf)
library(tidyverse)
library(gstat)

# SOURCE FUNCTIONS -----------------------------------

source("Functions/assignYearGroup().R")

# DATA FILES ------------------------------------------

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "../Data/Spatial_data/Boundaries_and_CRS/bng.prj")

bng <- st_crs("../Data/Spatial_data/Boundaries_and_CRS/bng.prj")$wkt

# Download UK boundary
UK_all <- geodata::gadm(country = "GBR", level = 2, 
                    path = "../Data/Spatial_data/Boundaries_and_CRS/UK")

# Climate data
climateDir <- "../Data/Spatial_data/Annual_climate_covar"

# Lakes file
lakesFile <- "../Data/Spatial_data/G2G/Lakes.tif"

# Woodland cover 1990 and 2015
LCM1990File <- "../Data/Spatial_data/LCM/1990/UK1990all_1km.tif"
LCM2015File <- "../Data/Spatial_data/LCM/2015/UK2015all_1km.tif"

# Connectivity test files
conn1990File <- "../Data/Spatial_data/Omniscape/radius160_resistance100/output1990/cum_currmap.tif"
conn2015File <- "../Data/Spatial_data/Omniscape/radius160_resistance100/output2015/cum_currmap.tif"

# SET PARAMETERS ------------------------------------

# List climate covariates
climateVar <- c("GDD5", "WMIN", "tasCV", "soilM", "RAIN")

# List all covariates
allVar <- c(climateVar, "coverBF", "coverCF", "connW")

# List range of years used for climate variables
range1 = c(1980, 2000)
range2 = c(2000, 2020)

### CODE --------------------------------------------

# READ COVARIATES -----------------------------------

# WOODLAND COVER COVARIATE

# Read in rasters
LCM1990 <- rast(LCM1990File) # 1990
LCM2015 <- rast(LCM2015File) # 2015

# Make broadleaf covariate
coverBF <- c(LCM1990["BF"], LCM2015["BF"])
names(coverBF) <- paste0(names(coverBF),"_", c("1990", "2015"))

# Make coniferous covariate
coverCF <- c(LCM1990["CF"], LCM2015["CF"])
names(coverCF) <- paste0(names(coverCF),"_", c("1990", "2015"))

# CONNECTIVITY COVARIATE

# Read in rasters and aggregate to 1km using median
conn1990 <- rast(conn1990File) %>%
  terra::aggregate(., fact = 40, fun = "median")
conn2015 <- rast(conn2015File) %>%
  terra::aggregate(., fact = 40, fun = "median")

# Make connectivity covariate
connW <- c(conn1990, conn2015)
names(connW) <-  paste0("conn_",c("1990", "2015"))

# Save aggregated spatRast for analysis
writeRaster(connW,
            "../Data/Spatial_data/Omniscape/omniConn_1km_forModel.tif",
            overwrite = TRUE)

# CLIMATE COVARIATES (AS MULTI-YEAR GROUP SPATRASTERS)

# Loop through climate variables
for (i in climateVar) {
  
  # Find annual file names for each climate covariate
  varFiles <- list.files(paste0(climateDir, "/", i),
                         full.names = TRUE,
                         pattern = "\\.tif$")
  
  # Extract year of files
  varFilesYear <- sub("\\.tif$", "", varFiles)  %>% # Drop ".tif"
    substr(., nchar(.) - 3 ,nchar(.)) # Extract last 4 'year' characters
  
  # Rasterise all annual rasters
  varR <- rast(varFiles)
  
  # Set names of spatRast to years
  names(varR) <- varFilesYear
  
  # Assign year groups (different ranges to species data)
  visitYears <- sapply(names(varR), function(x) {
    
    assignYearGroup(x,
                    yearGroupRange1 = range1,
                    yearGroupRange2 = range2) })
  
  # Find unique year group numbers (excluding NAs)
  visitYearGroups <-  unique(visitYears) %>%
    .[!is.na(.)]
  
  # For each multi-year group ...
  varGroupedR <- lapply( visitYearGroups, function(x) {
    
    # Subset all years to year group 'x', and take the mean
    varGroupedMean <- varR[[visitYears == x]] %>% 
      mean(.)
    
    # Change name to group name
    names(varGroupedMean) <- paste0("YearGroup", as.character(x))
    
    return(varGroupedMean)
    
  }) %>% rast(.) # Merge into same spatRaster
  
  # Assign to object
  assign(i, varGroupedR)
}

# PROCESS UK SPATIAL FILES ------------------------

# LAKES

# Create lake spatRast
# N.B. connW spatRast used as spatRast template
lakes_R <- rast(lakesFile) %>%
  crop(., connW) %>%
  extend(., connW)

# Create lakes spatVect object for later record filtering step
lakes <- as.polygons(lakes_R)

# UK BOUNDARY
# N.B. cannot use 'NAME_1' or 'NAME_2' columns to select countries as 
# 'geodata' package GADM issue, may be fixed in future.

# Remove Shetland and Orkney as extremely low tree cover
UK <- UK_all[UK_all$GID_1 != "GBR.3_1" | # Include all non-Scottish objects
               UK_all$GID_1 == "GBR.3_1" & # Also, if object is Scottish...
               UK_all$HASC_2 != "GB.OR",] # only include if not Orkney...

# Aggregate to remove boundaries
UK <- aggregate(UK)

# Convert to bng
UK <- project(UK, bng)

### Remove small islands

# Disaggregate
UK <- disagg(UK)

# Calculate area
UK$area_sqkm <- expanse(UK, unit = "km")

# Remove polygons with < 10km ^2 area
UK <- UK[UK$area_sqkm > 10]

# Aggregate back
UK <- aggregate(UK)

# UK SMOOTHED BOUNDARY
# Create smoothed boundary for use in modelling/mesh creation from UK

# Create spatRast for UK of any 1x1km cell where land cover >0
# N.B. connW spatRast used as spatRast template
coverUK_R <- rasterize(UK, connW,
                     touches = TRUE,
                     cover = TRUE)

# Filter UK spatRast to where land cover >50%
UK_R <- ifel(coverUK_R > 0.5,
             yes = 1,
             no = NA)

# Remove lake cells from UK_R
# N.B. This layer is used as prediction surface later on
UK_R <- mask(UK_R, lakes_R, inverse = TRUE)

# Create UK spatRast polygon from UK_R
UK_poly <- as.polygons(UK_R) %>% # Convert to polygon 
  st_as_sf(.) # Convert to sf object for smoothing

# Process UK_poly (smoothing, and convert back to spatVector)
# N.B. This is used to created mesh and functions as modelling boundary (domain)
smoothUK <- smoothr::smooth(UK_poly, method = "chaikin") %>%
  vect(.)

# INTERPOLATE  VARIABLES -------------------------------------

# Create binary raster of cells within and touching smoothUK,
# i.e. where we might need covariate values. 
# N.B. Holes need to be filled otherwise some records near lakes 
# or complex coastline have NA values
allUK_R <- st_as_sf(smoothUK) %>%#  Convert to sf object 
  smoothr::fill_holes(., threshold = Inf) %>% # Fill holes
  rasterize(., connW, # Rasterise using connW as template
            touches = TRUE)

# Loop through covariates
for (i in allVar) {

  # Get covariate spatRast
  iVar <- get(i)
  
  # Save layer names
  layerNames <- names(iVar)
  
  # Create data frame of covariate values and coordinates
  iVar_df <- as.points(iVar) %>%  
    data.frame(., crds(.))
  
  # For each layer name, i.e. each year...
  iVarInterpolated <- lapply(layerNames, function(j) {
    
    # Create a gstat formula for interpolation for layer j
    ijIntModel <- gstat(formula = get(j) ~ 1, # Only use location
                        locations = ~ x + y, # Location based on x,y
                        data = iVar_df, 
                        nmax = 5) # Only use nearest 5 points
    
    # Interpolate using gstat formula for layer j
    ijInterpolated <- interpolate(allUK_R,
                                  ijIntModel,
                                  na.rm = TRUE)[["var1.pred"]]
    
    }) %>% rast() # Join all layers together into one spatRast
  
  # Revert to original names
  names(iVarInterpolated) <- layerNames
  
  # Overwrite covariate object with interpolated object
  assign(i, iVarInterpolated)
  
}

# MASK COVARIATES TO UK DOMAIN -------------------------------

# Loop through covariates (and lakes spatRast)
for (i in allVar) {
  
  # Get raster object
  iVar <- get(i)
  
  # Mask to any 1x1km cell touching UK
  iVar <- terra::mask(iVar, allUK_R)
  
  # Align with connW spatRast template
  iVar <- crop(iVar, connW) %>%
    extend(., connW)
  
  # Overwrite object
  assign(i, iVar)
}

# CONVERT TO KM RESOLUTION -----------------------------------
# N.B. This is needed for INLA/inlabru

# VECTORS

# Reproject smoothUK to km for inlabru (needs changing back after!)
smoothUK <- project (smoothUK, gsub( "units=m", "units=km",
                                     st_crs(bng)$proj4string ))

# RASTERS

# Reproject covariates, and convert from sf to sp
for (i in c(allVar, "UK_R")) {

  # Reproject to km for inlabru
  assign(i, project (get(i), gsub( "units=m", "units=km",
                              st_crs(bng)$proj4string )))
}

# SCALE COVARIATES ------------------------------------------------------
# Need to scale fixed covariates for interpretation/interaction effects,
# and climate covariates to avoid very large/small numbers causing numerical issues (also comparable)

# Create scaling parameters data frame to populate
scalingParams <- data.frame(variable = allVar,
                            variableMean = rep(NA, length = length(allVar)),
                            variableSD = rep(NA, length = length(allVar)))

# Loop through each unscaled covariate
for (i in scalingParams$variable) {

  # Get spatRaster object
  iVar <- get(i)
  
  # Extract values from all iVar layers
  iVarValues <- iVar[] %>% # Extract values 
    as.matrix %>% # Convert into matrix
    as.vector %>% # Merge layer columns into single vector
    na.omit # Remove NAs
  
  # Calculate mean and SD using values from all layers
  iVarMean <- mean(iVarValues)
  iVarSD <- sd(iVarValues)
  
  # Save mean and SD to use later
  scalingParams[ scalingParams$variable == i, "variableMean"] <- iVarMean
  scalingParams[ scalingParams$variable == i, "variableSD"] <- iVarSD
  
  # Use parameters to scale the data
  iVar <- ( iVar - iVarMean ) / iVarSD

  # Overwrite object
  assign(paste0(i, "_scaled"), iVar)
}

# GROUP CLIMATE COVARIATES FOR SPLINE FITTING --------------------------

# Loop through climate variables
for (i in climateVar) {
  
  # Get climate SGDF i
  iVar <- paste0(i, "_scaled") %>% get(.)
  
  # Use the inla.group function to group covariate values
  # into 'n'
  iVar[] <- INLA::inla.group(iVar[],
                             n = 50,
                             method = "cut")
  
  # Assign to object
  assign(paste0(i, "_grp"), iVar)
}

# CREATE INTERACTION TERMS -------------------------------------------

# Broadleaf and coniferous interaction with connectivity
coverBF_connW <- coverBF_scaled * connW_scaled
coverCF_connW  <- coverCF_scaled * connW_scaled

# SAVE OBJECTS -------------------------------------------------------

# Create directories
dir.create("../Data/Spatial_data/DataForInlabru/spatVector",
           recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/Spatial_data/DataForInlabru/spatRaster",
           showWarnings = FALSE)

# Save scaling parameters
save(scalingParams,
     file = "../Data/Spatial_data/DataForInlabru/scalingParams.RData")

# Save UK vect objects
writeVector(UK,
            filename = "../Data/Spatial_data/DataForInlabru/spatVector/UK.shp",
            layer = "UK",
            overwrite = TRUE)
writeVector(smoothUK,
            filename = "../Data/Spatial_data/DataForInlabru/spatVector/smoothUK.shp",
            layer = "smoothUK",
            overwrite = TRUE)

# Save other spatRasts
for (i in c("UK_R", "coverBF", "coverCF", "connW",
          "GDD5", "WMIN", "tasCV", "soilM", "RAIN",
          "coverBF_scaled", "coverCF_scaled", "connW_scaled",
          "GDD5_scaled", "WMIN_scaled", "tasCV_scaled", "soilM_scaled", "RAIN_scaled",
          "GDD5_grp", "WMIN_grp", "tasCV_grp", "soilM_grp", "RAIN_grp",
          "coverBF_connW", "coverCF_connW")) {
  
  writeRaster(get(i),
              filename = paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                                i,
                                ".tif"),
                                overwrite = TRUE)
}
