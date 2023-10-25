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
library(raster)
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

# Read in rasters (need to aggregate to 1km)
conn1990 <- rast(conn1990File) %>%
  terra::aggregate(., fact = 40, sum, na.rm = TRUE) / (40^2)
conn2015 <- rast(conn2015File) %>%
  terra::aggregate(., fact = 40, sum, na.rm = TRUE) / (40^2)

# Make connectivity covariate
connW <- c(conn1990, conn2015)
names(connW) <-  paste0("conn_",c("1990", "2015"))

# Save aggregated spatRast for analysis
writeRaster(connW,
            "../Data/Spatial_data/Omniscape/omniConn_1km.tif",
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
    
    # Create a gstat formala for interpolation for layer j
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
for (i in c(allVar)) {
  
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

# CONVERT TO 'SP' AND KM RESOLUTION -----------------------------------
# N.B. This is needed for INLA/inlabru but creates warning messages

# VECTORS

# Convert smoothUK
smoothUK_sp <- as(st_as_sf(smoothUK), 'Spatial')
proj4string(smoothUK_sp) <- st_crs(bng)$proj4string

# Reproject smoothUK to km for inlabru (needs changing back after!)
smoothUK_sp <- spTransform(smoothUK_sp, gsub( "units=m", "units=km",
                                              st_crs(bng)$proj4string ))

# RASTERS

# Convert UK_R from spatRaster to raster
UK_R_sp <- stack(UK_R)

# Reproject to km for inlabru
extent(UK_R_sp) <- extent(c(xmin(UK_R_sp), xmax(UK_R_sp), ymin(UK_R_sp), ymax(UK_R_sp))/1000)
projection(UK_R_sp) <- gsub("units=m", "units=km", projection(UK_R_sp))

# SPATIAL GRID DATA FRAMES

# Reproject covariates, and convert from sf to sp
for (i in c(allVar)) {
  
  # Get spatRaster i
  iVar <- get(i)
  
  # Convert from spatRaster to raster
  iVar <- stack(iVar)
  
  # Reproject to km for inlabru
  extent(iVar) <- extent(c(xmin(iVar), xmax(iVar), ymin(iVar), ymax(iVar))/1000)
  projection(iVar) <- gsub("units=m", "units=km", projection(iVar))
  
  # Assign to new spatial grid data frameobject
    as(iVar, "SpatialGridDataFrame") %>%
      assign(paste0(i, "_unscaled_SGDF"), ., envir = .GlobalEnv)

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
  iVar <- paste0(i, "_unscaled_SGDF") %>%
    get(.)
  
  # Extract values from all iVar layers
  iVarValues <- iVar@data %>% # Extract SGDF values 
    as.matrix %>% # Convert into matrix
    as.vector %>% # Merge layer columns into single vector
    na.omit # Remove NAs
  
  # Calculate mean and SD using values from all layers
  iVarMean <- mean(iVarValues)
  iVarSD <- sd(iVarValues)
  
  # Save mean and SD to use later
  scalingParams[ scalingParams$variable == i, "variableMean"] <- iVarMean
  scalingParams[ scalingParams$variable == i, "variableSD"] <- iVarSD
  
  # Use parameters to scale the data (need to convert back to raster briefly)
  iVar <- ( stack(iVar) - iVarMean ) / iVarSD
  iVar <- as(iVar, "SpatialGridDataFrame")
  
  # Overwrite object
  assign(paste0(i, "_SGDF"), iVar)
}

# CREATE INTERACTION TERMS -------------------------------------------

# Broadleaf and coniferous interaction with connectivity
coverBF_connW_SGDF <- ( stack(coverBF_SGDF) * stack(connW_SGDF) ) %>%
  as(.,  "SpatialGridDataFrame")

coverCF_connW_SGDF  <- ( stack(coverCF_SGDF) * stack(connW_SGDF) ) %>%
  as(.,  "SpatialGridDataFrame")

# SAVE OBJECTS -------------------------------------------------------

# Save UK vect objects
writeVector(UK,
            filename = "../Data/Spatial_data/Boundaries_and_CRS/StudyArea(UK)/UK.shp",
            layer = "UK",
            overwrite = TRUE)
writeVector(smoothUK, 
            filename = "../Data/Spatial_data/Boundaries_and_CRS/StudyArea(UK)/smoothUK.shp",
            layer = "smoothUK",
            overwrite = TRUE)

# Save other spatial objects as .RData file
save(UK_R_sp, smoothUK_sp, scalingParams,
     coverBF_SGDF, coverCF_SGDF, connW_SGDF,
     coverBF_connW_SGDF, coverCF_connW_SGDF,
     GDD5_SGDF, WMIN_SGDF, tasCV_SGDF, soilM_SGDF, RAIN_SGDF,
     file = "../Data/Spatial_data/DataForInlabru.RData")
