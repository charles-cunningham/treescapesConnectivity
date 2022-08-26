# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Carabid inlabru spatio-temporal model test script
#
# Script Description: Test run
#
#

# TO DO LIST ----------------------------------------

# - Cover and connectivity to be tested separately and added in
# - Minimum winter temperature, change code to include end of previous year (i.e. Dec t-1, Jan t, Feb t)
# - Add min year parameter, probably 1990
# - Investigate climate smoothing terms
# - (Decompose out climate covariates into time and space)?
# - Change to another taxa group to test, not enough visits for carabids to test in National Forest!
# - 

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Install.packages
# install.packages("janitor")

# Load packages
library(INLA) 
library(inlabru)
library(terra)
library(GADMTools)
library(tidyverse)

# DATA FILES ------------------------------------------

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.prettywkt?download",
              destfile = "../Data/Spatial_data/bng.prj")

bng <- "../Data/Spatial_data/bng.prj"

# Read in country outline from GADM
UK <- gadm_sf_loadCountries( fileNames = "GBR", level = 1,
                             basefile = "../Data/Spatial_data/")$sf

# Set CRS
UK <- st_transform(UK, st_crs(bng))

# Load in Carabid test data
load("../Data/Species_data/Carabids_170316_Cleaned_Data.rdata")

# Cliamte data
climateDir <- "../Data/Spatial_data/CHESS/Current_climate_variables"

# Study area boundary polygon
testAreaFile <- "../Data/Spatial_data/NFC boundary/Forest boundary.shp"

# LOAD FUNCTIONS ------------------------------------

# Load grid coordinate to meters conversion function
source("Functions/gridCoords().R")

# SET PARAMETERS ------------------------------------

# Species to model
mySpecies <- "COL_6133"

# List climate covariates
climateVar <- c("GDD5", "MTCO", "tasCV", "soilM")

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 50 ### Set v coarse for testing, best value will probably be ~50

### CODE --------------------------------------------

# PROCESS GB SPATIAL FILES ------------------------

# Select countries

# Remove Northern Ireland (not enough data?)
GB <- UK[UK$NAME_1 != "Northern Ireland", "NAME_1"]

# Only Wales for testing?
#GB <- UK[UK$NAME_1 == "Wales", "NAME_1"]

# Write layer to shapefile
st_write(GB, "../Data/Spatial_data/GB.shp", append = FALSE)

# Read in as SpatVector
GB <- vect("../Data/Spatial_data/GB.shp", crs = bng)

# Disaggregate
GB <- disagg(GB)

# Calculate area
GB$area_sqkm <- expanse(GB, unit="km")

# Remove polygons with < 20km ^2 area
GB <- GB[GB$area_sqkm > 10]

# Aggregate back
GB <- aggregate(GB) #aggregate back

# CONVERT GRIDREF TO METERS ----------------------------

# Create x and y columns
taxa_data$X <- NA
taxa_data$Y <- NA

# Populate x and y columns using gridCoords function
# N.B. Function find the bottom, left-hand corner of the 1km grid square, 
# so need to add 500m
for (i in unique(taxa_data$TO_GRIDREF)) { # Loop through each row
  
  # Run gridCoords()
  coordData <- gridCoords(i, units = "m")
  
  # Add in x and y
  taxa_data[taxa_data$TO_GRIDREF == i, "X"] <- coordData$x + 500
  taxa_data[taxa_data$TO_GRIDREF == i, "Y"] <- coordData$y + 500
  
}

# PROCESS TO VISITS DATAFRAME -------------------------

# Add visit column (gridref plus date)
taxa_data$VISIT <- paste(taxa_data$TO_GRIDREF, taxa_data$TO_STARTDATE, sep="_")

# Add count of number of species recorded at each visit
# (Number of species reported on a visit, as an indicator of effort -
# most records are of list length 1 and probably not a comprehensive survey)
taxa_data <- taxa_data %>% 
  add_count(VISIT, name = "NUM_SP")

# Transform into visits data.frame
visitData <- taxa_data %>%
  add_column(PRESENCE = 1) %>% # Add a presence column
  pivot_wider(names_from = CONCEPT,
              values_from = PRESENCE) %>% # Transform to wide data frame (species names now columns)
  janitor::clean_names(case = "all_caps") %>% # All characters in data.frame to CAPS
  rename(PRESENCE = toupper(mySpecies)) %>% # Change the species of interest column name to 'Presence'
  dplyr::select(TO_GRIDREF,TO_STARTDATE,
                YEAR, VISIT, PRESENCE, 
                X, Y, NUM_SP) # Convert back to long data.frame

# Factorise visit species length
visitData$VISIT_LENGTH <- ifelse(visitData$NUM_SP == 1,"single",
                                 ifelse(visitData$NUM_SP %in% 2:3, "short", "long" ))

# Change non-detections visits to zeros
visitData$PRESENCE[is.na(visitData$PRESENCE)] <- 0

# Check number of presense and pseudo-absences
table(visitData$PRESENCE)

# CROP TO STUDY AREA ----------------------------------

# Read study area
studyArea <- vect(testAreaFile, crs = bng)

# Make visitData into sf object
visitDataSpatial <- visitData %>%
  st_as_sf(.,coords=c("X","Y"), crs = st_crs(bng))

# Plot
ggplot() +
  geom_point(data = subset(visitData, PRESENCE == 0), aes(x=X, y=Y),
             color = "white") +
  geom_point(data = subset(visitData, PRESENCE == 1), aes(x=X, y=Y),
             color = "black") +
  geom_sf(data = st_as_sf(GB), colour = "black", fill = NA)

# Some coordinates in the Irish sea so need to filter these

# Filter to study area (replace with GB/UK boundary for final run)
visitData <- visitData %>%
  filter(st_intersects(visitDataSpatial, st_as_sf(studyArea), sparse = FALSE)[,1]) 

# PROCESS COVARIATES -----------------------------------

# AGGREGATE ANNUAL DATA TO MULTI-YEAR GROUPS

# Year to nearest decade function (N.B. 1970 - 79, 1980 - 89 etc. not 1971-1980, 1981-1990)
floor_decade <- function(value){ return(value - value %% 10) }

# Find species data decade
visitData$decade <- floor_decade(visitData$YEAR)
visitData$iYear <- ((visitData$decade - min(visitData$decade)) / 10) + 1

# CREATE EFFORT COVARIATE

### Add day of year column ( will be included in model as covariate)
visitData$WOY <- visitData$TO_STARTDATE %>%
  strftime(., format = "%V") %>%
  as.numeric(.)

# Sort factor levels (need to convert to dummy variables)
visitData <- visitData %>%
  model.matrix(object = ~VISIT_LENGTH) %>%
  as.data.frame() %>%
  dplyr::select(-1) %>%
  cbind(visitData, .)

# READ CLIMATE COVARIATES AS MULTI-YEAR GROUP SPATRASTERS

# Loop through climate variables
for (i in climateVar) {
  
  # Find files names for each climate covariate
  varFiles <- list.files(paste0(climateDir, "/", i),
                         full.names = TRUE,
                         pattern = "\\.tif$")
  
  # Extract year of file
  varFilesYear <- sub("\\.tif$", "", varFiles)  %>% # Drop ".tif"
    substr(., nchar(.) - 3 ,nchar(.)) %>% # Extract last 4 'year' characters
    as.integer(.) # Convert to number
  
  # Convert to decade as per species data above
  varFilesYearGroups <- floor_decade(varFilesYear)
  
  # Find multi-year groups with visit records (if 0 for some study areas, need to remove)
  visitYears <- varFilesYearGroups[varFilesYearGroups %in% visitData$decade] %>% unique(.)
  
  # For each multi-year group (decade)...
  varGroupedR <- lapply( visitYears, function(x) {
    
    # Rasterize associated files
    varRast <- varFiles[varFilesYearGroups == x] %>%
      rast(.)
    
    # Take the mean of the multi-year group
    varGroupedMean <- mean(varRast)
    
    # Change name to group name
    names(varGroupedMean) <- as.character(x)
    
    return(varGroupedMean)
    
  }) %>% rast(.) # Merge into same spatRaster
  
  # Assign to object
  assign(paste0(i, "_R"), varGroupedR)
  
}

### SCALE COVARIATES?

#...


# CONVERT TO 'SP' AND KM RESOLUTION FOR INLABRU (CUE WARNING MESSAGES) -----

# Load raster package for 'sp' stages
library(raster)

# Convert visitData
visitData_sp <- visitData
coordinates(visitData_sp) <- c("X", "Y")
proj4string(visitData_sp) <- st_crs(bng)$proj4string

# Convert GB
studyArea_sp <- as(st_as_sf(studyArea), 'Spatial')
proj4string(studyArea_sp) <- st_crs(bng)$proj4string

# Reproject to km for inlabru (needs changing back after!)
studyArea_sp <- spTransform(studyArea_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))
visitData_sp <- spTransform(visitData_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# Reproject all covariates, and convert from sf to sp (currently only climate)
for (i in climateVar) {
  
  # Get climate raster
  climVar <- paste0(i, "_R") %>% get(.)
  
  # Convert from spatRaster to raster
  climVar <- stack(climVar)
  
  # Reproject to km for inlabru
  climVar <- projectRaster(climVar, crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))
  
  # Assign to new object
  assign(paste0(i, "_R"), climVar)
  
}

# INTERPOLATE COVARIATE NA VALUES AND CREATE COVARIATE FUNCTIONS --------

# INTERPOLATE

# Loop through climate variables
for (i in climateVar) {

  # Get climate raster
  climR <- paste0(i, "_R") %>% get(.)
  
  # Extract covariate data for visit locations
  visitCovData <- extract(climR, visitData_sp)
  
  # Extract NA values (i.e. which visit locations are NA for any time step)
  visitCovNA <- visitData_sp[rowSums(is.na(visitCovData)) > 0,]
  
  # Check if any NA values, if 0...
  if (length(visitCovNA) == 0) {
    
    # ...just convert climR to SpatialGridDataFrame and assign
    as(climR, "SpatialGridDataFrame") %>%
      assign(paste0(i, "_SGDF"), ., envir = .GlobalEnv)
    
  } else { # ...otherwise, need to interpolate NA values
    
    # Convert coordinates of NA values to spatial points object
    NApoints <- visitCovNA %>% 
      coordinates(.) %>% # Extract coordinates
      unique(.)%>% # Only unique values
      SpatialPoints(., proj4string = crs(visitData_sp)) # Convert to SpatialPoints object
    
    # Interpolate around cells (buffer = 2 means nearest 13 cells are included (mean))
    visitCovInterpolate <- extract(x = climR, y = NApoints, buffer = 2, fun = mean)
    
    # Create SPDF of interpolated cells
    interpolateSPDF <- SpatialPointsDataFrame(coords = coordinates(NApoints),
                                              data = data.frame(visitCovInterpolate),
                                              proj4string = crs(visitData_sp))
    
    # Create duplicate raster to add interpolated cells onto
    climR_interpolated <- climR
    climR_interpolated[] <- NA
    
    # Update raster with interpolated cells, each row seperately
    for(j in 1:nlayers(climR)) {
      
      climR_interpolated[[j]] <- rasterize(interpolateSPDF[j], climR[[j]], names(interpolateSPDF[j]), update=TRUE)
    }
    
    # Convert to SpatialGridDataFrame and assign
    as(climR_interpolated, "SpatialGridDataFrame") %>%
      assign(paste0(i, "_SGDF"), ., envir = .GlobalEnv)
  }
}

# CREATE COVARIATE FUNCTIONS

# fun(x,y)

f.soilS <- function(x, y) {
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(soilM_SGDF[1]))
  proj4string(spp) <- fm_sp_get_crs(soilM_SGDF[1])
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, soilM_SGDF[1])
 
  return(v[,1])
}

# fun(x,y,time)

# This is tricky to construct,N.B.
# Finn Lindgren says "For inlabru to use the function properly, it needs to be vectorised, so that
# fun(x,y,year) allows x, y, and year to be vectors, and returns a vector v of elements v_i, where
# v_i = fun(x[i], y[i], year[i]) that is, one value for each row of data.frame(x, y, year)"
f.soilST <- function(x, y, iYear) {
  
  # turn coordinates into SpatialPoints object:
  # with the appropriate coordinate reference system (CRS)
  spp <- SpatialPointsDataFrame(coords = data.frame(x = x, y = y), 
                       data = data.frame(iYear = iYear),
                       proj4string = fm_sp_get_crs(soilM_SGDF))
  
  # Extract elevation values at spp coords, from our SpatialGridDataFrame
  v <- over(spp, soilM_SGDF)
  
  v$iYear <- iYear
  
  v$iYearCov <- apply(v, 1, FUN = function(x) { x[ x["iYear" ]] })
  
  return(v$iYearCov)
}

# Test
f.soilST(coordinates(visitData_sp)[,"X"],coordinates(visitData_sp)[,"Y"], visitData_sp$iYear )

# CREATE MESH -------------------------------------------------------

# Max edge is range/8 as a rule of thumb
maxEdge = estimated_range/10 # estimated_range/8 for Wales

mesh <- inla.mesh.2d(boundary = studyArea_sp,
                     loc = visitData_sp,
                     max.edge = c(1,2) * maxEdge, #  "max.edge = c(1,4) * maxEdge," for Wales, but crashes for National forest
                     cutoff = maxEdge/2,
                     offset = c(1,2) * maxEdge, # "max.edge = c(1,4) * maxEdge," for Wales, but crashes for National forest
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

ggplot() + 
  gg(mesh) + 
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE == 0)), color = "white") +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE == 1)), color = "black") +
  coord_fixed() +
  geom_sf(data = st_as_sf(studyArea_sp), col = "black", fill = NA)

# FIT SPATIO-TEMPORAL MODEL ---------------------------------

# Create indices
iYear <- visitData_sp$iYear
nYear <- length(unique(visitData_sp$iYear))

# Define spatial pattern
mySpace <- inla.spde2.pcmatern(
  mesh,
  alpha = 2,
  prior.range = c(10 * maxEdge, 0.5),   
  prior.sigma = c(1, 0.5))

# Set components
inlabruCmp  <-  PRESENCE ~
  WOY(main = WOY , model = "ar1") + 
  soilST(main = f.soilST(X, Y, iYear), model = "linear") + # change to smooth term?
  VISIT_LENGTHshort + VISIT_LENGTHsingle +
  spatial(main = coordinates,
          group = iYear,
          ngroup = nYear,
          model = mySpace,
          control.group=list(model="ar1")) +
  Intercept(1)
  
# Fit model
inlabruST2 <- bru(components = inlabruCmp,
                 family = "binomial",
                 data = visitData_sp,
                 options=list(verbose=TRUE))
#Model summary
summary(inlabruST2)

# PREDICT -----------------------------------------

# Create spatial data frame to predict over
ppxl <- pixels(mesh, mask = studyArea_sp)
ppxlAll <- cprod(ppxl, data.frame(iYear = seq_len(nYear)))

# Predict using spatio-temporal model
# ( N.B. exlcuding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
inlabruPred <- predict(inlabruST2, 
                       ppxlAll, 
                       ~ data.frame(iYear = iYear,
                                    lambda =  plogis(Intercept + 
                                                       spatial +
                                                       soilST + #Change to model used
                                                       max(inlabruST2$summary.random$WOY$mean)))) # max value for WOY to predict over
# Plot
ggplot() + 
  gg(inlabruPred, aes(fill=mean)) +
  scale_fill_viridis_c(option="magma",direction=1) +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE==1)), color = "blue") +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE==0)), color = "red") +
  facet_wrap(~iYear) +
  coord_fixed() +
  geom_sf(data = st_as_sf(studyArea_sp), col = "black", fill = NA)

# Add plots to look at posteriors and fixed effects
