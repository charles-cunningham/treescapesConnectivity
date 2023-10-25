# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Inlabru spatio-temporal model script
#
# Script Description: ...
#
# R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

# TODO ----------------------------------------

# - Set covariate priors
# - (Decompose out climate covariates into time and space)?
# - Add in northern ireland
# - Final term/interaction formulation decision
# - Final connectivity metric decision
# - Evaluation plots

# # INSTALL PACKAGES -----------------------------------
# #Run this code once
# 
# # Install INLA
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/stable"),
#                  dep=TRUE)
# 
# # Install inlabru
# remotes::install_github("inlabru-org/inlabru", ref = "stable")
# 
# # Install sparta
# remotes::install_github("BiologicalRecordsCentre/sparta")
# 
# # Needed to run on HPC
# inla.binary.install()
# 
# # Other packages as required from CRAN, i.e install.packages()

# LOAD PACKAGES -----------------------------------

library(INLA) 
library(inlabru)
library(terra) # terra' must be loaded before 'sf' due to 'geodata' dependency bug
library(sf)
library(tidyverse)
library(raster)

# When running on cluster implement PARDISO
inla.setOption(pardiso.license = "Treescapes/pardiso.lic")
inla.pardiso.check()

# SOURCE FUNCTIONS ----------------------------

source("Treescapes/ST_SDMs/Functions/assignYearGroup().R")

# SET PARAMETERS ---------------------------------------

# Set number of threads
inla.setOption("num.threads" = 2)

# Each iteration is one species, specified array in job script
args <- commandArgs(trailingOnly = TRUE)
batchN <- as.numeric(args[1])

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 50 # if 150 works try 50 

# DATA FILES ------------------------------------------

### SPATIAL DATA

# Processed spatial files
load("Treescapes/ST_SDMs/Data/DataForInlabru.RData")

# Study area boundary polygon
# National Forest: "../Data/Spatial_data/Boundaries_and_CRS/NFC/Forest_boundary.shp"
# Mersey Forest:  "../Data/Spatial_data/Boundaries_and_CRS/EnglandCommunityForests/Englands_Community_Forests_Nov2021.shp"
testAreaFile <- "Treescapes/ST_SDMs/Data/Forest_boundary.shp"

### Download BNG WKT string
# N.B. Individual filename needed for each task to prevent
# different tasks tring to read/write at same time 

# Specify temporary file to download 'bng' CRS wkt to
tempFile <- paste0("Treescapes/ST_SDMs/Data/",batchN, "bng.prj")

# Download file
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = tempFile)

# Assign wkt string
bng <- sf::st_crs(tempFile)$wkt

# Remove temporary file
unlink(tempFile)

### SPECIES DATA

# Load in test data (using butterflies for now)
taxaData <- read.csv("Treescapes/ST_SDMs/Data/ccunningham_butterflies.csv")

### CODE --------------------------------------------

# PROCESS TO VISITS DATAFRAME -------------------------

### FILTER RAW DATA

# Assign raw data
taxaData <- rawData

# How many species/aggregates to model within taxa
# unique(taxaData$taxon) %>%
#   length(.)
# = 65

# Set species for this job
iSpecies <- unique(taxaData$taxon)[batchN]

# Set filename-friendly species name
iSpeciesTidy <- gsub(" ", "_", iSpecies) %>%
  gsub("-", ".", .)

# Filter to GB for now
taxaData <- taxaData[taxaData$srs == "OSGB", ]

# Remove 5x5km records that have crept into some data requests
taxaData <- taxaData[grep(pattern = '\\NE$|\\SE$|\\SW$|\\SE$',
                          x = taxaData$monad,
                          invert = TRUE) , ]

# Remove multi-day records as not usable in this analysis
taxaData <- taxaData[grep(pattern = '-',
                          x = taxaData$date,
                          invert = TRUE) , ]

# Remove any duplicate rows (shouldn't be any but just in case!)
# N.B. Each row is a record and there can be multiple species records per visit.
taxaData <- unique(taxaData)

### ADD DATA COLUMNS

# Add visit column (concatenate gridref plus date)
taxaData$visit <- paste(taxaData$monad, taxaData$date, sep = "_")

# Add iSpecies presence column
taxaData$presence <- if_else(taxaData$taxon == iSpecies, 1, 0)

# Add count of number of records at each visit
# (Number of records reported on a visit, as an indicator of effort -
# most records are of list length 1 and probably not a comprehensive survey)
taxaData <- add_count(taxaData, visit, name = "numRecords")

# Factorise visit species length
taxaData$visitLength <- ifelse( taxaData$numRecords == 1,
                                "Single",
                                ifelse(taxaData$numRecords %in% 2:3,
                                       "Short", "Long"))

# Select useful columns
taxaData <- dplyr::select(taxaData,
                          monad, date, visit, presence, numRecords, visitLength)

### REDUCE ROW COUNT TO UNIQUE VISITS

# Filter all records to a single value for each visit
visitData <- taxaData %>%
  # For each visit ...
  group_by(visit) %>%
  # Extract max presence value
  slice(which.max(presence)) %>%
  # Finally, ungroup
  ungroup()

# Check number of presence and pseudo-absences
table(visitData$presence)

# CONVERT GRIDREF TO METERS ----------------------------
# N.B. this may also be a more convenient option: https://github.com/colinharrower/BRCmap

# Convert grid reference to Latitude and Longitude
visitLatLong <- sparta::gr2gps_latlon(visitData$monad, 
                                      precision = rep(1000,length(visitData$monad )),  
                                      projection = rep("OSGB", length(visitData$monad )),
                                      centre = TRUE)

# Convert to Easting and Northing
visitXY <-  sparta::gps_latlon2gr(latitude = visitLatLong[,"LATITUDE"],
                                  longitude = visitLatLong[,"LONGITUDE"],
                                  out_projection = "OSGB",
                                  return_type = "en" )

# Round tp nearest 500m as slight (~1m) rounding errors in previous two conversions
visitXY$EASTING <- round(visitXY$EASTING / 500) * 500
visitXY$NORTHING <- round(visitXY$NORTHING / 500) * 500

# Change column names to x and y
colnames(visitXY) <- c("x", "y")

# Join together
visitData <- bind_cols(visitData, visitXY)

# Remove data frames no longer needed 
rm(visitLatLong, visitXY)

# CROP SPECIES RECORDS TO STUDY AREA --------------------------

# Read study area (DROP for actual run)
GB <- testAreaFile %>%
  vect(., crs = bng)

# Make visitData into sf object
visitDataSpatial <- visitData %>%
  st_as_sf(., coords=c("x","y"), crs = bng)

# Filter to study area
visitData <- visitData %>%
  filter(st_intersects(visitDataSpatial, st_as_sf(GB), sparse = FALSE)[,1]) 

# PROCESS COVARIATES -----------------------------------

# AGGREGATE ANNUAL DATA TO MULTI-YEAR GROUPS

# Add year column
visitData$year <- visitData$date %>%
  as.POSIXct(., format = "%d/%m/%Y") %>%
  format(., format="%Y") %>% 
  as.numeric(.)

### Add week of year column ( will be included in model as covariate)
# N.B. Week of the year as decimal number (00--53) using Monday as 
# the first day of week (and typically with the first Monday 
# of the year as day 1 of week 1). The UK convention.
visitData$week <- visitData$date %>%
  as.POSIXct(., format = "%d/%m/%Y") %>%
  strftime(., format = "%W") %>%
  as.numeric(.)

# Find species year group
visitData$iYear <- sapply(visitData$year, assignYearGroup)

# Drop years outside year groups
visitData <- visitData[!is.na(visitData$iYear),]

# CREATE EFFORT COVARIATE

# Convert factor levels to dummy variables
visitData <- visitData %>%
  model.matrix(object = ~visitLength) %>%
  as.data.frame() %>%
  dplyr::select(-1) %>%
  cbind(visitData, .)

# CONVERT TO 'SP' AND KM RESOLUTION FOR INLABRU (CUE WARNING MESSAGES) -----

# Convert visitData
visitData_sp <- visitData
coordinates(visitData_sp) <- c("x", "y")
proj4string(visitData_sp) <- st_crs(bng)$proj4string

# Reproject to km for inlabru (needs changing back after!)
visitData_sp <- spTransform(visitData_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# Convert GB
GB_sp <- as(st_as_sf(GB), 'Spatial')
proj4string(GB_sp) <- st_crs(bng)$proj4string

# Reproject vectors to km for inlabru (needs changing back after!)
GB_sp <- spTransform(GB_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# CREATE MESH -------------------------------------------------------

# Max edge is range/8 as a rule of thumb (range/3 to range/10)
maxEdge <- estimated_range/5 # estimated_range/8 for Wales, this is a very sensitive parameter : estimated_range/4 works

# Find record locations to build mesh from
recordCoords <- coordinates(visitData_sp) %>% 
  unique(.)

mesh <- inla.mesh.2d(boundary = GB_sp,   # smoothGB_sp,
                     loc = recordCoords,
                     max.edge = c(1,5) * maxEdge, #  "max.edge = c(1,4) * maxEdge," for Wales, but crashes for National forest
                     offset = c(1,2) * maxEdge, # "max.edge = c(1,4) * maxEdge," for Wales, but crashes for National forest
                     cutoff = maxEdge/10,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# FIT SPATIO-TEMPORAL MODEL ---------------------------------

# Create indices
iYear <- visitData_sp$iYear
nYear <- length(unique(visitData_sp$iYear))

# Define spatial SPDE priors
mySpace <- inla.spde2.pcmatern(
  mesh,
  alpha = 2,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(2, 0.5))

# Fixed effect priors
C.F. <- list( mean = 0,
              prec = 1 ) # Precision for all fixed effects except intercept

# RW priors
pcprec <- list(prior = "pcprec",
               param = c(1, 0.01))

# Set components
inlabruCmp  <-  presence ~
  soilST(main = soilM_SGDF_grp,
         main_selector = "iYear",
         model = "rw2",
         scale.model = TRUE,
         hyper = list(prec = list(param = pcprec$param))) +
  wminST(main = WMIN_SGDF_grp,
         main_selector = "iYear",
         model = "rw2",
         scale.model = TRUE,
         hyper = list(prec = list(param = pcprec$param))) +
  cvST(main = tasCV_SGDF_grp,
       main_selector = "iYear",
       model = "rw2",
       scale.model = TRUE,
       hyper = list(prec = list(param = pcprec$param))) +
  gddST(main = GDD5_SGDF_grp,
        main_selector = "iYear",
        model = "rw2",
        scale.model = TRUE,
        hyper = list(prec = list(param = pcprec$param))) +
  coverBF(main = coverBF_SGDF,
          main_selector = "iYear",
          model = "linear") +
  coverCF(main = coverCF_SGDF,
          main_selector = "iYear",
          model = "linear") +
  connectivity(main = connW_SGDF,
               main_selector = "iYear",
               model = "linear") +
  BFconnINT(main = coverBF_connW_SGDF,
            main_selector = "iYear",
            model = "linear") +
  CFconnINT(main = coverCF_connW_SGDF,
            main_selector = "iYear",
            model = "linear") +
  visitLengthShort(main = visitLengthShort,
                   main_selector = "iYear",
                   model = "linear") +
  visitLengthSingle(main = visitLengthSingle,
                    main_selector = "iYear",
                    model = "linear") +
  week(main = week,
       main_selector = "iYear",
       model = "rw1",
       scale.model = TRUE,
       cyclic = TRUE,
       hyper = list(prec = list(param = pcprec$param))) +
  spatial(main = coordinates,
          group = iYear,
          ngroup = nYear,
          model = mySpace,
          control.group = list(model = "ar1")) +
  Intercept(1)

# Fit model
inlabruST <- bru(components = inlabruCmp,
                 family = "binomial",
                 control.family = list(link = "cloglog"),
                 data = visitData_sp, 
                 options=list(control.fixed = C.F.,
                              control.inla= list(int.strategy='eb'),
                              verbose = TRUE))

# Output model summary
summary(inlabruST)

# PREDICT -----------------------------------------

### Create grid prediction pixels
ppxl <- as(mask(GB_R_sp, GB_sp), "SpatialPixelsDataFrame")

# Create multi-time period prediction pixels
ppxlAll <- cprod(ppxl, data.frame(iYear = seq_len(nYear)))

# Predict using spatio-temporal model
# ( N.B. exlcuding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
inlabruPred <- predict(inlabruST, 
                       ppxlAll, 
                       ~ data.frame(iYear = iYear,
                                    lambda =  plogis( spatial +
                                                        soilST +
                                                        wminST +
                                                        cvST +
                                                        gddST +
                                                        coverBF +
                                                        coverCF +
                                                        connectivity +
                                                        BFconnINT +
                                                        CFconnINT +
                                                        max(inlabruST$summary.random$week$mean) + # max value for WOY to predict over
                                                        Intercept ))) # Long visit is inferred

# SAVE -----------------------------------------

# Assign objects to species specific name (model fit, summary, prediction)
assign(paste0(iSpeciesTidy, "_STmodel"), inlabruST)
assign(paste0(iSpeciesTidy, "_STmodelSummary"), summary(inlabruST))
assign(paste0(iSpeciesTidy, "_STpred"), inlabruPred)

# Save objects(model fit, summary, prediction)
save(list = paste0(iSpeciesTidy,
                   "_STmodel"),
     file = paste0(
       "Treescapes/ST_SDMs/Output/LocalTest/",
       iSpeciesTidy,
       ".fit.Rdata"))

save(list = paste0(iSpeciesTidy,
                   "_STmodelSummary"),
     file = paste0(
       "Treescapes/ST_SDMs/Output/LocalTest/",
       iSpeciesTidy,
       ".summary.Rdata"))

save(list = paste0(iSpeciesTidy,
                   "_STpred"),
     file = paste0(
       "Treescapes/ST_SDMs/Output/LocalTest/",
       iSpeciesTidy,
       ".pred.RData"))

##########################################################################################

# ### N.B. COVARIATE FUNCTIONS INFO
# 
# ## (i) X-Y (SPACE) COVARIATE FUNCTION
# 
# # Option 1: Function to map out single-layer spatial grid data frame (soilM_SGDF[1])
# 
# f.soilS <- function(x, y) {
#   
#   # turn coordinates into SpatialPoints object with the appropriate coordinate reference system (CRS)
#   spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(soilM_SGDF[1]))
#   
#   # Extract values at spp coords
#   v <- over(spp, soilM_SGDF[1])
#   
#   return(v[,1])
# }
# 
# # And component detailed as:
# 
# ... + soil(f.soilS(x, y), model = "linear")
# 
# 
# # Option 2:
# 
# # Enter directly as component using inbuilt spatial grid data frame recognition
# 
# ... + soil(soilM_SGDF[1], model = "linear")
# 
# ## (i) X-Y-TIME (SPACE-TIME) COVARIATE FUNCTION
# 
# # Option 1: Function to map out multi-layer spatial grid data frame (each layer is time step)(soilM_SGDF)
# 
# # This is tricky to construct, N.B. Finn Lindgren says
# # "For inlabru to use the function properly, it needs to be vectorised, so that
# # fun(x,y,year) allows x, y, and year to be vectors, and returns a vector v of elements v_i, where
# # v_i = fun(x[i], y[i], year[i]) that is, one value for each row of data.frame(x, y, year)"
# 
# f.soilST <- function(x, y, iYear) {
# 
#   # turn coordinates into SpatialPoints object with the appropriate coordinate reference system (CRS)
#   spp <- SpatialPointsDataFrame(coords = data.frame(x = x, y = y),
#                                 data = data.frame(iYear = iYear),
#                                 proj4string = fm_sp_get_crs(soilM_SGDF))
# 
#   # Extract values at spp coords, from our SpatialGridDataFrame
#   v <- over(spp, soilM_SGDF)
# 
#   v$iYear <- iYear
# 
#   v$iYearCov <- apply(v, 1, FUN = function(z) { z[ z["iYear" ]] })
# 
#   return(v$iYearCov)
# }

# # And component detailed as:
# 
# ... + soilST(main = f.soilST(x, y, iYear), model = "linear")
# 
# 
# # Option 2:
# 
# # Enter directly as component using inbuilt spatial grid data frame recognition
# 
# ... + soilST(main = soilM_SGDF, main_selector = "iYear", model = "linear")
