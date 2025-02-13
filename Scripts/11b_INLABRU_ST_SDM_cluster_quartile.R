# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Inlabru spatio-temporal model script
#
# Script Description: Same as script 11a, but subsets sample area to broadleaf cover quantiles
# and omit all cover from analysis 
#
# R version 4.3.2 (2023-10-31) -- "Eye Holes"
# Copyright (C) 2023 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

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
# install.packages("inlabru")
#
# # Install BRCmap (more complex to install as NAMESPACE issue - 
# rgdal and rgeos included but not on CRAN)
# Work around (hopefully this will be fixed soon): Download locally, extract from .zip, 
# remove rgdal and rgeos from NAMESPACE file, re -zip, then run code below -
# devtools::install_local("Packages/BRCmap-master.zip")
# install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
# 
# # Needed to run on HPC
# inla.binary.install()
# 
# # Other packages as required from CRAN, i.e install.packages()

# LOAD PACKAGES -----------------------------------

library(INLA) 
library(inlabru)
library(sf)
library(terra)
library(stars)
library(tidyverse)
library(BRCmap)
library(ggplot2)
library(grid)
library(gridExtra)

# SOURCE FUNCTIONS ----------------------------

source("Treescapes/ST_SDMs/Functions/assignYearGroup().R")

# SET PARAMETERS ---------------------------------------

# Organise raw data taxa groups
dataBC <- c( "Butterflies", "Moths" )
dataBRC <- c( "Bryophytes", "Carabids", "Caddisflies", "Centipedes",
              "Ephemeroptera", "Gelechiidae", "Hoverflies", "Ladybirds",
              "Lichen", "Molluscs", "Odonata", "Orthoptera",
              "Shieldbugs", "Soldierflies", "Spiders")

# List range of years used for species records
range1 = c(1990,2000)
range2 = c(2015,2025)

# Set number of threads
numThreads = 4

# Set batch number/species (as.numeric(args[1])) and taxa group (arg[2]), specified array in job script
args <- commandArgs(trailingOnly = TRUE)
batchN <- as.numeric(args[1])
taxaGroup <- args[2]

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 50 

# Set bitmap type for ggplot2::ggsave to work on cluster
options(bitmapType='cairo')

# SET UP INLA -----------------------------------------

# Set number of threads
inla.setOption("num.threads" = numThreads)
# ...and print out
paste(numThreads, "threads used") %>%
  print

# DATA FILES ------------------------------------------

### SPATIAL DATA

# Scaling parameters
load("Treescapes/ST_SDMs/Data/DataForInlabru/scalingParams.RData")

# SpatRasters
for (i in list.files("Treescapes/ST_SDMs/Data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("Treescapes/ST_SDMs/Data/DataForInlabru/spatRaster/",
                     i))) 
}

# SpatVectors
for (i in list.files("Treescapes/ST_SDMs/Data/DataForInlabru/spatVector",
                     pattern =  "\\.shp$")) {
  
  assign(gsub(".shp", "", i),
         vect(paste0("Treescapes/ST_SDMs/Data/DataForInlabru/spatVector/",
                     i)))
}

### Download BNG WKT string
# N.B. Individual filename needed for each task to prevent
# different tasks tring to read/write at same time 

# Specify temporary file to download 'bng' CRS wkt to
tempFile <- paste0("Treescapes/ST_SDMs/Data/",taxaGroup, batchN, "bng.prj")

# Download file
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = tempFile)

# Assign wkt string
bng <- sf::st_crs(tempFile)$wkt

# Remove temporary file
unlink(tempFile)

### SPECIES DATA

# LOAD IN SPECIES DATA (according to taxa group)

# Bryophytes (only taxa with separate NI and GB files, so need to process)
if(taxaGroup == "Bryophytes") {
  
  # Read in GB and NI files separately
  rawDataGB <- read.csv("Treescapes/ST_SDMs/Data/Raw_data/Bryophytes/Treescapes_Bryophytes_220822.csv")
  rawDataNI <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Bryophytes/Bryophytes_2022.rds")
  
  # Rename GB monad and species columns
  rawDataGB <- rename(rawDataGB, 
                      monad = Monad,
                      recommended_name = AGG_NAME)
  
  # Process GB date columns
  rawDataGB <- mutate(rawDataGB,
                      startdate = paste0(YEAR1, "-", MONTH1, "-", DAY1),
                      enddate = paste0(YEAR1, "-", MONTH1, "-", DAY1))
  
  # Join NI and GB together
  rawDataUK <- bind_rows(rawDataGB, # Bind all rows together
                         rawDataNI) %>%
    dplyr:: select(all_of(intersect(names(rawDataGB), # Select only shared columns
                                    names(rawDataNI)))) } 

# Butterflies
if(taxaGroup == "Butterflies") {
  rawDataUK <- read.csv("Treescapes/ST_SDMs/Data/Raw_data/Butterflies/ccunningham_butterflies.csv") }

# Carabids
if (taxaGroup == "Carabids") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Carabids/Ground_Beetles_2022.rds")}

# Caddisflies
if (taxaGroup == "Caddisflies") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Caddisflies/Caddisflies_2022_Connected Treescapes.rds") }

# Centipedes
if (taxaGroup == "Centipedes") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Centipedes/Centipedes_2022.rds") } 

# Ephemeroptera
if (taxaGroup == "Ephemeroptera") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Ephemeroptera/Ephemeroptera_2022.rds") }

# Gelechiidae
if (taxaGroup == "Gelechiidae") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Gelechiidae/Gelechiidae_2022.rds") }

# Hoverflies
if (taxaGroup == "Hoverflies") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Hoverflies/Hoverflies_2022.rds") }

# Ladybirds
if (taxaGroup == "Ladybirds") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Ladybirds/Ladybirds_2022.rds") }

# Lichen
if (taxaGroup == "Lichen") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Lichen/Lichen_2022.rds") }

# Molluscs (need to remove marine species)
if (taxaGroup == "Molluscs") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Molluscs/Molluscs_2022_Connected Treescapes.rds")
  
  # Filter marine species out of dataset
  rawDataUK  <- rawDataUK %>%
    filter(taxon_group_one != "marine molluscs"  | is.na(taxon_group_one)) %>%
    filter(taxon_group_two != "marine molluscs"  | is.na(taxon_group_two)) }

# Moths
if (taxaGroup == "Moths") {
  rawDataUK <- read.csv("Treescapes/ST_SDMs/Data/Raw_data/Moths/ccunningham_moths.csv") }

# Odonata
if (taxaGroup == "Odonata") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Odonata/Odonata_2022.rds") }

# Orthoptera
if (taxaGroup == "Orthoptera") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Orthoptera/Orthoptera_2022.rds") }

# Shieldbugs
if (taxaGroup == "Shieldbugs") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Shieldbugs/Shieldbugs_2022.rds") }

# Soldierflies
if (taxaGroup == "Soldierflies") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Soldierflies/Soldierflies_2022.rds") }

# Spiders
if (taxaGroup == "Spiders") {
  rawDataUK <- readRDS("Treescapes/ST_SDMs/Data/Raw_data/Spiders/Spiders_2022.rds") }

# FILTER AND STANDARDISE DATA FRAMES --------------------------

### ASSIGN RAW DATA

# Assign raw data
taxaData <- rawDataUK

### STANDARDISE COLUMN NAMES BETWEEN TAXA GROUPS

# Process species name (DEPENDS ON TAXA GROUP)
if(taxaGroup %in% dataBRC) {
  
  # Remove subgenus "(subgenus)" from a taxon name
  taxaData$recommended_name <- trimws(gsub("\\([[:alpha:]]*\\)[[:space:]]", 
                                           "", 
                                           taxaData$recommended_name))
  # Standardise species name column
  taxaData <- rename(taxaData, 
                     taxon = recommended_name)
}

### FILTER RAW DATA
# N.B. Each row is a record and there can be multiple species records per visit.

# Filter by monad

# Remove records with no monad data
taxaData <- taxaData[ !is.na(taxaData$monad), ]

# Remove 5x5km records that may have crept into some data requests
taxaData <- taxaData[grep(pattern = '\\NE$|\\SE$|\\SW$|\\SE$',
                          x = taxaData$monad,
                          invert = TRUE) , ]

# Filter by date

# Format date (DEPENDS ON TAXA GROUP)
if(taxaGroup %in% dataBC) {
  
  # Remove multi-day records as not usable in this analysis
  taxaData <- taxaData[grep(pattern = '-',
                            x = taxaData$date,
                            invert = TRUE) ,]
  # Reformat date
  taxaData$date <- as.POSIXct(taxaData$date, format = "%d/%m/%Y")
}

if(taxaGroup %in% dataBRC) {
  
  # Remove multi-day records as not usable in this analysis
  taxaData <- taxaData[taxaData$startdate == taxaData$enddate, ]
  
  # Create new standardised "date" column or, if column already exists, overwrite
  taxaData <- mutate(taxaData, date = startdate)
  
  # Reformat date
  taxaData$date <- as.POSIXct(taxaData$date, format = "%Y-%m-%d")
}

## Add year column
taxaData$year <- taxaData$date %>%
  format(., format="%Y") %>% 
  as.numeric(.)

# Aggregate annual data to multi-year groups
taxaData$iYear <- sapply(taxaData$year, function(x) {
  
  assignYearGroup(x,
                  yearGroupRange1 = range1,
                  yearGroupRange2 = range2) })

# Drop years outside year groups
taxaData <- taxaData[!is.na(taxaData$iYear),]

# Select useful columns
taxaData <- dplyr::select(taxaData, taxon, monad, date, year, iYear)

# PROCESS TO VISITS DATAFRAME ---------------------------------

### GET SPECIES NAME

# Set species for this model
iSpecies <- unique(taxaData$taxon)[batchN]

# How many species/aggregates to model within taxa?
unique(taxaData$taxon) %>%
  length %>%
  paste("Total taxon to model within group:", . ) %>%
  print

### N.B. Total taxon to model within different groups: --------
# Butterflies = 64
# Bryophytes = 969 #
# Caddisflies = 186
# Carabids = 371 
# Centipedes = 45
# Ephemeroptera = 54
# Gelechiidae = 151
# Hoverflies = 262
# Ladybirds = 53
# Lichen = 2047 #
# Molluscs = 208
# Moths = 864
# Odonata = 51
# Orthoptera = 65
# Shieldbugs = 75
# Soldierflies = 154
# Spiders = 674

# Print for run log
paste0("Species: ", iSpecies,
      ", batchN: ", batchN) %>%
  print

# Set filename-friendly species name
iSpeciesTidy <- gsub(" ", "_", iSpecies) %>%
  gsub("-", ".", .)

### ADD DATA COLUMNS

# Add visit column (concatenate gridref plus date)
taxaData$visit <- paste(taxaData$monad, taxaData$date, sep = "_")

# Add iSpecies presence column
taxaData$presence <- if_else(taxaData$taxon == iSpecies, 1, 0)

# Drop taxon column
taxaData <- dplyr::select(taxaData, -taxon)

# Add count of number of records at each visit
# (Number of records reported on a visit, as an indicator of effort -
# most records are of list length 1 and probably not a comprehensive survey)
taxaData <- add_count(taxaData, visit, name = "numRecords")

# Factorise visit species length
taxaData$visitLength <- ifelse( taxaData$numRecords == 1,
                                "Single",
                                ifelse(taxaData$numRecords %in% 2:3,
                                       "Short", "Long"))

### REDUCE ROW COUNT TO UNIQUE VISITS

# Filter all records to a single value for each visit
visitData <- taxaData %>%
  # For each visit ...
  group_by(visit) %>%
  # Extract max presence value (i.e. if any records are iSpecies then 1, else 0)
  slice(which.max(presence)) %>%
  # Finally, ungroup
  ungroup()

# Check number of presence and pseudo-absences
table(visitData$presence)

# CONVERT GRIDREF TO EASTING/NORTHING --------------------------

# Add coordinate reference system ("OSGB" or "OSNI")
visitData$crs <- gr_det_country(visitData$monad)

# Remove data if crs not "OSNI" or "OSGB" (i.e drop channel islands and NAs)
visitData <- visitData[visitData$crs %in% c("OSNI", "OSGB"),]

# Convert OS monads to OSGB easting and northing
# N.B. Irish OSNI monads are on a different projection so are not multiples of 500
visitXY <- OSgrid2GB_EN(visitData$monad,
                        centre = TRUE,
                        gr_prec = 1000)

# Add to visitData
visitData <- bind_cols(visitData, visitXY)

# Create column to add recentred OSNI easting/northing to
visitData$OSGB_NORTHING <- visitData$OSGB_EASTING <- NA

# Recentre OSNI easting/northing to centre of OSGB grid
visitData <- visitData %>%
  mutate(OSGB_EASTING = case_when(crs == "OSGB" ~ 
                                    EASTING,
                                  crs == "OSNI" ~ 
                                    floor( EASTING/1000) * 1000 + 500)) %>%
  mutate(OSGB_NORTHING = case_when(crs == "OSGB" ~
                                     NORTHING,
                                   crs == "OSNI" ~ 
                                     floor(NORTHING/1000) * 1000 + 500))

# Rename columns and remove ones not needed further
visitData <- visitData %>%
  rename(x = OSGB_EASTING, # Rename coordinates to x and y
         y = OSGB_NORTHING) %>%
  dplyr::select(-c(crs, # Remove columns not needed further
                   EASTING, NORTHING))

# CROP SPECIES RECORDS TO STUDY AREA --------------------------

# Make sure in km
smoothUK <- project (smoothUK, 
                     gsub( "units=m", "units=km",
                           st_crs(bng)$proj4string ))

# Convert to terra vector object
visitDataSpatial <- vect(visitData, geom=c("x", "y"), crs=bng, keepgeom=TRUE)

# Reproject to km from m (needed to avoid inla numeric issues)
visitDataSpatial <- project (visitDataSpatial,
                             gsub( "units=m", "units=km",
                                   st_crs(bng)$proj4string ))

# Mask records to within vector
visitDataSpatial <- mask(visitDataSpatial, smoothUK)

# PROCESS COVARIATES -----------------------------------

# CREATE WEEK COVARIATE

### Add week of year column ( will be included in model as covariate)
# N.B. Week of the year as decimal number (01--53) as defined in ISO 8601.
# If the week (starting on Monday) containing 1 January has four or more 
# days in the new year, then it is considered week 1. Otherwise, 
# it is the last week of the previous year, and the next week is week 1.
visitDataSpatial$week <- visitDataSpatial$date %>%
  strftime(., format = "%V") %>%
  as.numeric(.)

# CREATE EFFORT COVARIATE

# Convert factor levels to dummy variables
visitDataSpatial <- visitDataSpatial %>%
  model.matrix(object = ~visitLength) %>%
  as.data.frame() %>%
  dplyr::select(-1) %>%
  cbind(visitDataSpatial, .)

# FILTER DATA TO COVER BRACKETS -------------------------------------

# Set quantile 'brackets' in percentage [broadleaf cover only]
# (exclude landscapes with no broadleaf woodland)
coverBins <- cbind(c(0, 25),
                   c(25, 50),
                   c(50, 75),
                   c(75, 100))

# Set bracket names
colnames(coverBins) <- c("0-25", "25-50", "50-75", "75-100")

# Print quartiles (proportion) 
coverBF$BF_1990[] %>%
  .[!is.na(.)] %>%
  .[. > 0] %>%
  quantile(., probs = seq(0,1,0.25)) 

### SET LOOP HERE (runs until end of script)
for (coverBin in colnames(coverBins)) {

# Set quantile limits for bracket i using broadleaf cover [1990 only]
quantLimits <- coverBF$BF_1990[] %>%
  .[!is.na(.)] %>%
  .[. > 0] %>%
  quantile(., probs = coverBins[, coverBin]/100)

# Mask broadleaf cover (1990) to only cells that fall within bracket
quantCells_R <- ifel(coverBF$BF_1990 >= quantLimits[1] & 
                     coverBF$BF_1990 < quantLimits[2],
                   1,
                   NA)

# Create polygon version of quantCells
quantCells_V <- as.polygons(quantCells_R)

# Mask visit data to cells that fall within broadleaf cover bracket
visitDataQuant <- mask(visitDataSpatial, quantCells_V)

# Mask other data to cells that fall within broadelaf cover bracket
for (i in c( "GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp",
             "coverBF_scaled", "coverCF_scaled", "connW_scaled")) {
  
  assign(paste0(i, "Quant"), mask(get(i), quantCells_R))
  
}

# EXTRACT COVARIATES FOR EFFECTS PLOT -------------------------------
# Need to extract spatial covariates over species records for plots

### Create data frame of covariate values at visit locations

# Create a base data frame with iYear, presence and week (non-spatial) to build up from
covarValues <- dplyr::select(as.data.frame(visitDataQuant), iYear, presence, week)

# Loop through spatial (random) variables
for (i in c( "GDD5_grpQuant", "WMIN_grpQuant", "tasCV_grpQuant", "RAIN_grpQuant", "soilM_grpQuant",
             "coverBF_scaledQuant", "coverCF_scaledQuant", "connW_scaledQuant")) {

  # Get covariate i spatRaster (each layer is for time period iYear)
  cov_R <- get(i)
  
  # Extract values for all iYear layers and add to beginning of data frame using species records
  covarValues <- terra::extract(cov_R, visitDataQuant, ID = FALSE) %>%
    cbind(covarValues) # Add existing data frame onto end
  
  # Calculate correct covariate for iYear for each value
  # (i.e. record is only for one time period - we need that one)
  cov_iYear <- apply(covarValues, 1, FUN = function(x) { x[ x[ "iYear" ]] })
  
  # Add newly created correct iYear variable values to end of data frame, removing suffix
  covarValues[, gsub("_.*","", i)] <- cov_iYear
  
  # Drop now redundant columns from first columns
  covarValues <- covarValues[, -c(1,2)]
  
}

### Process climate covariate values for rug plot

# Remove cover and climate variables and convert to long format
climCovarValues <-  dplyr::select(covarValues,
                                  -c(coverBF, coverCF, connW)) %>%
  gather(randomEff, value, 3:NCOL(.))

# Count number of each quantile value for each variable for presence/absence separately
climCovarValues <- climCovarValues %>%
  group_by(randomEff, presence, value) %>%
  mutate(count = n()) %>%
  ungroup %>%
  distinct

# TIDY MEMORY BEFORE MODEL RUN --------------------------------------

# Remove data frames no longer needed 
rm(visitXY, rawDataUK, taxaData, visitData)

# Garbage clean
gc()

# CREATE MESH -------------------------------------------------------

# Max edge is as a rule of thumb (range/3 to range/10)
maxEdge <- estimated_range/5

# Find record locations to build mesh from
recordCoords <- crds(visitDataQuant) %>% 
  unique(.)

# Create mesh (convert boundary to sp object as leads to best convergence)
mesh <- inla.mesh.2d(boundary = st_as_sf(quantCells_V) %>% as("Spatial"),
                     loc = recordCoords,
                     max.edge = c(1,5) * maxEdge,
                     offset = c(1,2) * maxEdge, 
                     cutoff = maxEdge/2,
                     min.angle = 26,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# FIT SPATIO-TEMPORAL MODEL ---------------------------------

# Create indices
iYear <- visitDataQuant$iYear
nYear <- length(unique(iYear))

# Define spatial SPDE priors
mySpace <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(1, 0.5))

# Priors for fixed effects
fixedHyper <- list( mean = 0,
                    prec = 1 ) # Precision for all fixed effects except intercept

# Priors for random effects
randomHyper <- list(theta = list(prior="pc.prec",
                               param=c(0.5, 0.01)))
ar1Hyper <- list(rho = list(prior="pc.prec",
                            param=c(0.5, 0.01)))

# Set components
inlabruCmp  <-  presence ~ 0 + Intercept(1) +
  
  GDD5(main = GDD5_grpQuant,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  WMIN(main = WMIN_grpQuant,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  tasCV(main = tasCV_grpQuant,
        main_layer = iYear,
        model = "rw2",
        scale.model = TRUE,
        hyper = randomHyper) +
  RAIN(main = RAIN_grpQuant,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  soilM(main = soilM_grpQuant,
        main_layer = iYear,
        model = "rw2",
        scale.model = TRUE,
        hyper = randomHyper) +
  connectivity(main = connW_scaledQuant,
               main_layer = iYear,
               model = "linear") +
  visitLengthShort(main = visitLengthShort,
                   model = "linear") +
  visitLengthSingle(main = visitLengthSingle,
                    model = "linear") +
  week(main = week,
       model = "rw2",
       cyclic = TRUE,
       hyper = randomHyper) +
  spaceTime(main = geometry,
            group = iYear,
            ngroup = nYear,
            model = mySpace,
            control.group = list(model = "ar1",
                                 hyper = ar1Hyper))

# Fit model
model <- bru(components = inlabruCmp,
             family = "binomial",
             control.family = list(link = "cloglog"),
             data = st_as_sf(visitDataQuant),
             options=list(control.fixed = fixedHyper,
                          control.inla= list(int.strategy='eb'),
                          control.compute = list(waic = TRUE, dic = FALSE, cpo = TRUE),
                          verbose = TRUE))

# Assign model summary object and output
modelSummary <- summary(model); modelSummary

# PREDICT -----------------------------------------

# Create grid prediction pixels
ppxl <- mask(quantCells_R, smoothUK) %>%
  crop(.,smoothUK ) %>%
  as.points %>%
  st_as_sf

# Create multi-time period prediction pixels
ppxlAll <- fm_cprod(ppxl, data.frame( iYear = seq_len(nYear)))

# Predict using spatio-temporal model
# ( N.B. excluding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
modelPred <- predict(model, 
                     ppxlAll, 
                     ~ data.frame(
                                  iYear = iYear,
                                  lambda =  1 - exp( -exp( spaceTime + # cloglog back transform
                                                             soilM +
                                                             WMIN +
                                                             tasCV +
                                                             GDD5 +
                                                             RAIN +
                                                             connectivity +
                                                             # Max value for week to predict over (removed later)
                                                             max(model$summary.random$week$mean) + 
                                                             Intercept ))),
                     exclude = c("week")) 

# MODEL EVALUATION --------------------------

# SET PARAMETERS

# Labels
randomEffLabels <- c('GDD5' = "Growing degree days", 'RAIN' = "Annual precipiation",
                     'soilM' = "Soil moisture", 'tasCV' = "Temperature seasonality",
                     'week' = "Week of year" , 'WMIN' = "Winter minimum temperature")
timeLabels <- data.frame(label = c("1990 - 2000", "2015 -"), 
                         iYear = c("1", "2"))
linearEffLabels <- c('connectivity' = "Connectivity",
                     'visitLengthSingle' = "Single record visit",
                     'visitLengthShort' = "Short visit (2-3 records)")

# Template raster for converting from sf to terra raster objects
template_R <- st_as_stars(UK_R)
template_R[[1]][1:ncell(template_R)] <- NA

# CALCULATE LOGCPO

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

# NON-SPATIAL RANDOM EFFECTS PLOT

### Create effects data frame

# Extract random effects from model, and exclude spatial
randomEff_df <- model$summary.random
randomEff_df["spaceTime"] <- NULL

# Add name of random effect to each dataframe in list
randomEff_df <- imap(randomEff_df, ~mutate(.x, randomEff = .y))

# Unlist, then rename and select quantile columns
randomEff_df <- do.call(rbind, randomEff_df)%>%
  rename("q0.025" = "0.025quant",
         "q0.5" = "0.5quant",
         "q0.975" = "0.975quant") %>%
  dplyr::select(ID, q0.025, q0.5, q0.975, randomEff)

### Back scale non-spatial random covariate values

# Apply 'unscaling' function to every row of effects data frame
randomEff_df$unscaledID <- sapply(1:NROW(randomEff_df), function(x) {
  
  # If week, just use ID as not scaled
  if (randomEff_df$randomEff[x] == "week") {
    
    return(randomEff_df$ID[x])
    
  } else { # Otherwise, 'unscale'!
    
    # Extract covariate mean and sd for scaling function
    randomEffMean <-
      scalingParams[scalingParams$variable == randomEff_df$randomEff[x],
                    "variableMean"]
    randomEffSD <-
      scalingParams[scalingParams$variable == randomEff_df$randomEff[x],
                    "variableSD"]
    
    # Unscale using xSCALED = (x - xbar)/sd --> x = (xSCALED * sd) + xbar principle
    unscaledID <- (randomEff_df$ID[x] * randomEffSD) + randomEffMean
    
    return(unscaledID) # Return unscaled value
    
  }})

# Apply 'unscaling' function to every row of record locations
climCovarValues$unscaledValue <- sapply(1:NROW(climCovarValues), function(x) {
  
  # If week, just use value as not scaled
  if (climCovarValues$randomEff[x] == "week") {
    
    return(climCovarValues$value[x])
    
  } else { # Otherwise, 'unscale'!
    
    # Extract covariate mean and sd for scaling function
    randomEffMean <- scalingParams[scalingParams$variable == climCovarValues$randomEff[x],
                                   "variableMean"]
    randomEffSD <- scalingParams[scalingParams$variable == climCovarValues$randomEff[x],
                                 "variableSD"]
    
    # Unscale using xSCALED = (x - xbar)/sd --> x = (xSCALED * sd) + xbar principle
    unscaledValue <- (climCovarValues$value[x] * randomEffSD) + randomEffMean
    
    return(unscaledValue)
    
  }})

### Plot

randomEffPlot <- ggplot(randomEff_df) +
  
  # Random effect size
  geom_line(aes(x = unscaledID, y = q0.5)) +
  geom_line(aes(x = unscaledID, y = q0.025), lty = 2, alpha = .5) +
  geom_line(aes(x = unscaledID, y = q0.975), lty = 2, alpha = .5) +
  
  # Record covariate values for absence (bottom) and presence (top)
  geom_rug(data = subset(climCovarValues, presence == 0),
           aes(x = unscaledValue, alpha = count), linewidth = 1,
           sides = "b", show.legend = FALSE, inherit.aes = FALSE) +
  geom_rug(data = subset(climCovarValues, presence == 1),
           aes(x = unscaledValue, alpha = count), linewidth = 1,
           sides = "t", show.legend = FALSE, inherit.aes = FALSE) +
  
  # Thematics
  facet_wrap(~ randomEff, scale = 'free_x', labeller = as_labeller(randomEffLabels))  + 
  ggtitle("Non-linear random effects") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10)) +
  xlab("") +
  ylab("Probability of occurrence")

# FIXED EFFECTS PLOT

# Loop through covariates and extract estimates
for (i in names(linearEffLabels)) {
  
  # For covariate i, extract effect size
  effectSize <- modelSummary$inla$fixed[i,] %>% 
    t %>% # Transpose
    data.frame 
  
  # Add covariate
  effectSize$Covariate <- i
  
  # If first covariate
  if( i == names(linearEffLabels)[1]) {
    
    # Create a new data frame
    effectSizeAll <- effectSize
    
  }  else {
    
    # Join data frames together
    effectSizeAll <- rbind(effectSizeAll, effectSize) 
    
  }
}

# Plot fixed effects
fixedEffPlot <- ggplot(effectSizeAll, 
                       aes(y = X0.5quant, x = Covariate,
                           ymin = X0.025quant, ymax=X0.975quant, 
                           col = Covariate, fill = Covariate)) + 
  #specify position here
  geom_linerange(linewidth=4, colour = "lightblue") +
  ggtitle("Linear effects") +
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=2, shape=21, colour="white", fill = "black", stroke = 0.1) +
  scale_x_discrete(name="",
                   limits = rev(names(linearEffLabels)),
                   labels = as_labeller(linearEffLabels)) +
  scale_y_continuous(name="Effect size") +
  coord_flip() +
  theme_minimal() + 
  guides(colour = "none") +
  theme(axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, vjust = -0.5))

# SPDE PARAMETER POSTERIORS

# Extract range and variance of space-time SPDE, and plot
range.plot <- plot( spde.posterior(model, "spaceTime", what = "range")) +
  ggtitle("SPDE range") +
  theme(plot.title = element_text(hjust = 0.5))
var.plot <- plot(spde.posterior(model, "spaceTime", what = "log.variance")) +
  ggtitle("SPDE log variance") +
  theme(plot.title = element_text(hjust = 0.5))

# MATERN CORRELATION AND COVARIANCE

corplot <- plot(spde.posterior(model, "spaceTime", what = "matern.correlation")) +
  ggtitle("Matern correlation") +
  theme(plot.title = element_text(hjust = 0.5))
covplot <- plot(spde.posterior(model, "spaceTime", what = "matern.covariance")) +
  ggtitle("Matern covariance") +
  theme(plot.title = element_text(hjust = 0.5))

# MEDIAN AND SD PREDICTION PLOTS

### Plot median

# Convert to medians to spatRast for saving/plotting
median_R <- c(st_rasterize(modelPred[modelPred$iYear == "1", "median"],
                           template = template_R,
                           options = c("a_nodata = NA")),
            st_rasterize(modelPred[modelPred$iYear == "2", "median"],
                         template = template_R,
                         options = c("a_nodata = NA"))) %>%
  rast

# Add names, i.e. iYear
names(median_R) <- c("1", "2")

# Convert to data frame for plotting
median_df <- as.data.frame(median_R, xy = TRUE) %>% # Convert to data frame
  gather(iYear, median, c("1","2" )) # Convert to long format

# Plot posterior median
predMedian <- ggplot(data = median_df) +
  ggtitle("Median posterior occupancy") +
  coord_fixed() +
  geom_tile(aes(x=x, y=y, fill = median, colour = median)) +
  scale_fill_distiller(palette = "BuGn",
                       direction = 1,
                       limits = c(0,1),
                       guide = guide_colourbar(title = "Occupancy\nprobability")) +
  scale_colour_distiller(palette = "BuGn",
                         direction = 1,
                         limits = c(0,1),
                         guide = "none") +
  facet_wrap(~ iYear, labeller = as_labeller(timeLabels)) +
  geom_text(data = timeLabels,
            mapping = aes(x = 500, y = 925, label = label)) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        strip.text.x = element_blank()) +
  geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black")


### Plot posterior sd

# Convert to medians to spatRast for saving/plotting
sd_R <- c(st_rasterize(modelPred[modelPred$iYear == "1", "sd"],
                       template = template_R,
                       options = c("a_nodata = NA")),
            st_rasterize(modelPred[modelPred$iYear == "2", "sd"],
                         template = template_R,
                         options = c("a_nodata = NA"))) %>%
  rast

# Add names, i.e. iYear
names(sd_R) <- c("1", "2")

# Convert to data frame for plotting
sd_df <- as.data.frame(sd_R, xy = TRUE) %>% # Convert to data frame
  gather(iYear, sd, c("1","2" )) # Convert to long format

predSD <- ggplot(data = sd_df) +
  ggtitle("Posterior standard deviation") +
  geom_tile(aes(x=x, y=y, fill = sd, colour = sd)) +
  scale_fill_distiller(palette = "BuGn",
                       direction = 1,
                       guide = guide_colourbar(title = "Standard\ndeviation")) +
  scale_colour_distiller(palette = "BuGn",
                         direction = 1,
                         guide = "none") +
  facet_wrap(~ iYear, labeller = as_labeller(timeLabels)) +
  geom_text(data = timeLabels,
            mapping = aes(x = 500, y = 925, label = label)) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        strip.text.x = element_blank()) +
  coord_fixed() +
  geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black")

# POSTERIOR MEDIAN OF SPACE-TIME RANDOM FIELD

# Predict spatio-temporal field only at link scale
spaceTimePred <- predict(model, 
                         ppxlAll, 
                         ~ data.frame(iYear = iYear,
                                      effectSize =  spaceTime ),
                         include = c("spaceTime"))

# Convert to medians to spatRast for saving/plotting
spaceTime_R <- c(st_rasterize(spaceTimePred[spaceTimePred$iYear == "1", "median"],
                              template = template_R,
                              options = c("a_nodata = NA")),
                 st_rasterize(spaceTimePred[spaceTimePred$iYear == "2", "median"],
                              template = template_R,
                              options = c("a_nodata = NA"))) %>%
  rast

# Add names, i.e. iYear
names(spaceTime_R) <- c("1", "2")

# Convert to data frame for plotting
spaceTime_df <- as.data.frame(spaceTime_R, xy = TRUE) %>% # Convert to data frame
  gather(iYear, median, c("1","2" )) # Convert to long format

# Plot
spaceTimePlot <- ggplot(data = spaceTime_df) +
  ggtitle("Spatio-temporal field")  +
  geom_tile(aes(x=x, y=y, fill = median, colour = median)) +
  scale_fill_distiller(palette = "RdYlBu",
                       direction = 1,
                       guide = guide_colourbar(title = "SPDE\nposterior\nmedian"),
                       limits = c(-1,1) * max(abs(spaceTime_df$median))) +
  scale_colour_distiller(palette = "RdYlBu",
                         direction = 1,
                         guide = "none",
                         limits = c(-1,1) * max(abs(spaceTime_df$median))) +
  facet_wrap(~ iYear, labeller = as_labeller(timeLabels)) +
  geom_text(data = timeLabels,
            mapping = aes(x = 500, y = 925, label = label)) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        strip.text.x = element_blank()) +
  coord_fixed() +
  geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black", inherit.aes = FALSE)

### JOINT EVALUATION PLOT

evalPlot <- arrangeGrob(predMedian, 
                        predSD ,
                        fixedEffPlot ,randomEffPlot,
                        spaceTimePlot,
                        nrow = 4, ncol = 2,
                        layout_matrix = rbind(c(1, 1),
                                              c(2, 2),
                                              c(3, 4),
                                              c(5, 5)),
                                              top = grid::textGrob(paste0(iSpecies, ", 
                                                                          logCPO = ", logCPO),
                                                                   gp = grid::gpar(fontsize=20)))

# Compose matern plot
spdeAndMaternPlot <- arrangeGrob(range.plot, covplot,
                                 var.plot, corplot, ncol = 2)

# SAVE -----------------------------------------

# Define species directory
iSpeciesDir <- paste0("Treescapes/ST_SDMs/Output/Supplementary_quantile_analysis/", 
                      coverBin,
                      "_percent_cover_quartile/",
                      taxaGroup,
                      "/", 
                      iSpeciesTidy)

# Create species directory
if (!dir.exists(iSpeciesDir)) {
  
  dir.create(iSpeciesDir,
             recursive = TRUE)
}

### Save objects

# Model
save(model,
     file = paste0(iSpeciesDir,
                   "/modelFit.RData"))
# Model summary
save(modelSummary,
     file = paste0(iSpeciesDir,
                   "/modelSummary.RData"))

# Model prediction object
save(modelPred,
     file = paste0(iSpeciesDir,
                   "/modelPred.RData"))

# Posterior median relative occurrence probability prediction
writeRaster(median_R,
            file = paste0(iSpeciesDir,
                          "/medianPred.tif"),
            overwrite= TRUE)

# Evaluation plot
ggsave(paste0(iSpeciesDir,
              "/effectsPlot_range_", estimated_range,
              "_logCPO_", logCPO, ".png"),
       evalPlot,
       width = 4000, height = 8000, 
       units = "px", dpi = 400,
       limitsize = FALSE)

# SPDE parameter posterior (range and variance) and
# matern correlation and covariance plot
ggsave(paste0(iSpeciesDir,
              "/MatCorCovPlot_range_", estimated_range,
              "_logCPO_", logCPO, ".png"),
       spdeAndMaternPlot,
       width = 6000, height = 3000,
       units = "px", dpi = 400,
       limitsize = FALSE)

}
