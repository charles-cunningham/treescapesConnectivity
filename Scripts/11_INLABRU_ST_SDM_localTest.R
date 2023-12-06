# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Inlabru spatio-temporal model script
#
# Script Description: Takes raw recording data set and processes, and then fits a 
# SDM using inlabru. This version fits a spatio-temporal model with relative
# occurrence probability as the response using a visits data structure.
#
# R version 4.3.2 (2023-10-31 ucrt) -- "Eye Holes"
# Copyright (C) 2023 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)

# TODO ----------------------------------------

# - Consider rw2 or inla.mesh.1d for climate covariates

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to R: (C: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# INSTALL PACKAGES -----------------------------------
# #Run this code once
#
# # Install INLA
# install.packages("INLA", repos=c(getOption("repos"),
#                                  INLA="https://inla.r-inla-download.org/R/stable"),
#                  dep=TRUE)
# 
# # Install inlabru
# install.packages("inlabru")
#
# # Need terra >1.7.60 to work as fixed extract issue (issue 1332) (not on CRAN at time of writing)
# install.packages("terra", repos = "https://rspatial.r-universe.dev")
#
# # Install BRCmap (more complex to install as NAMESPACE issue - 
# rgdal and rgeos included but not on CRAN)
# Work around (hopefully this will be fixed soon): Download locally, extract from .zip, 
# remove rgdal and rgeos from NAMESPACE file, re -zip, then run code below -
# devtools::install_local("R:/rsrch/cb751/lab/Charles/R/PackageLibrary/BRCmap-master.zip")
# 
# Needed to run on HPC
# inla.binary.install()
#
# # Other packages as required from CRAN, i.e install.packages()

# LOAD PACKAGES ---------------------------------------

library(terra)
library(sf)
library(stars)
library(INLA) 
library(inlabru)
library(tidyverse)
library(BRCmap)
library(ggplot2)
library(grid)
library(gridExtra)

# When running on cluster implement PARDISO
# inla.setOption(pardiso.license = "Treescapes/pardiso.lic")
# inla.pardiso.check()

# SOURCE FUNCTIONS ----------------------------

source("Functions/assignYearGroup().R")

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
#numThreads = 4

# Set batch number/species (arg[1]) and taxa group (arg[2]), specified array in job script
#args <- commandArgs(trailingOnly = TRUE)
batchN <- 11
taxaGroup <- "Butterflies"

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 20 

# Set bitmap type for ggplot2::ggsave to work on cluster
#options(bitmapType = 'cairo')

# DATA FILES ------------------------------------------

### SPATIAL DATA

# Scaling parameters
load("../Data/Spatial_data/DataForInlabru/scalingParams.RData")

# Study area boundary polygon
# National Forest: "../Data/Spatial_data/Boundaries_and_CRS/NFC/Forest_boundary.shp"
# Mersey Forest:  "../Data/Spatial_data/Boundaries_and_CRS/EnglandCommunityForests/Englands_Community_Forests_Nov2021.shp"
smoothUK <- terra::vect("../Data/Spatial_data/Boundaries_and_CRS/NFC/Forest_boundary.shp")

# SpatRasters
for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
 
  assign(gsub(".tif", "", i),
         rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                                i)))
}

### Download BNG WKT string
# N.B. Individual filename needed for each task to prevent
# different tasks tring to read/write at same time 

# Specify temporary file to download 'bng' CRS wkt to
tempFile <- paste0("..//",batchN, "bng.prj")

# Download file
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = tempFile)

# Assign wkt string
bng <- sf::st_crs(tempFile)$wkt

# Remove temporary file
unlink(tempFile)

### SPECIES DATA

# Load in test data (butterflies for now)
rawDataUK <- read.csv("../Data/Species_data/Raw_data/Butterflies/ccunningham_butterflies.csv")

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
# Moths = 864
# Bryophytes = 969
# Carabids = 371 
# Caddisflies = 186
# Centipedes = 45
# Ephemeroptera = 54
# Gelechiidae = 151
# Hoverflies = 262
# Ladybirds = 53
# Lichen = 2047
# Molluscs = 280
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
# N.B. Week of the year as decimal number (00--53) using Monday as 
# the first day of week (and typically with the first Monday 
# of the year as day 1 of week 1). The UK convention.
visitDataSpatial$week <- visitDataSpatial$date %>%
  strftime(., format = "%W") %>%
  as.numeric(.)

# CREATE EFFORT COVARIATE

# Convert factor levels to dummy variables
visitDataSpatial <- visitDataSpatial %>%
  model.matrix(object = ~visitLength) %>%
  as.data.frame() %>%
  dplyr::select(-1) %>%
  cbind(visitDataSpatial, .)

# EXTRACT COVARIATES FOR EFFECTS PLOT -------------------------------
# Need to extract spatial covariates over species records for plots

### Create data frame of covariate values at visit locations

# Create a base data frame with iYear, presence and week (non-spatial) to build up from
covarValues <- dplyr::select(as.data.frame(visitDataSpatial), iYear, presence, week)

# Loop through spatial (random) variables
for (i in c( "GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp",
             "coverBF_scaled", "coverCF_scaled", "connW_scaled")) {

  # Get covariate i spatRaster (each layer is for time period iYear)
  cov_R <- get(i)
  
  # Extract values for all iYear layers and add to beginning of data frame using species records
  covarValues <- terra::extract(cov_R, visitDataSpatial, ID = FALSE) %>%
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

# CREATE MESH -------------------------------------------------------

# Max edge is as a rule of thumb (range/3 to range/10)
maxEdge <- estimated_range/5

# Find record locations to build mesh from
recordCoords <- crds(visitDataSpatial) %>% 
  unique(.)

# Create mesh
mesh <- inla.mesh.2d(boundary = st_as_sf(smoothUK),
                     loc = recordCoords,
                     max.edge = c(1,5) * maxEdge,
                     offset = c(1,2) * maxEdge, 
                     cutoff = maxEdge/5,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# FIT SPATIO-TEMPORAL MODEL ---------------------------------

# Create indices
iYear <- visitDataSpatial$iYear
nYear <- length(unique(iYear))

# Define spatial SPDE priors
mySpace <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(1, 0.5))

# Priors for fixed and random effects
fixedHyper <- list( mean = 0,
                    prec = 1 ) # Precision for all fixed effects except intercept
randomHyper <- list(theta = list(prior="pc.prec",
                                 param=c(0.5, 0.01)))
ar1Hyper <- list(rho = list(prior="pc.prec",
                            param=c(0.5, 0.01)))

# Set components
inlabruCmp  <-  presence ~ 0 + Intercept(1) +
  soilM(main = soilM_grp,
        main_layer = iYear,
        model = "rw2",
        scale.model = TRUE,
        hyper = randomHyper) +
  WMIN(main = WMIN_grp,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  tasCV(main = tasCV_grp,
        main_layer = iYear,
        model = "rw2",
        scale.model = TRUE,
        hyper = randomHyper) +
  GDD5(main = GDD5_grp,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  RAIN(main = RAIN_grp,
       main_layer = iYear,
       model = "rw2",
       scale.model = TRUE,
       hyper = randomHyper) +
  coverBF(main = coverBF_scaled,
          main_layer = iYear,
          model = "linear") +
  coverCF(main = coverCF_scaled,
          main_layer = iYear,
          model = "linear") +
  connectivity(main = connW_scaled,
               main_layer = iYear,
               model = "linear") +
  BFconnINT(main = coverBF_connW,
            main_layer = iYear,
            model = "linear") +
  CFconnINT(main = coverCF_connW,
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
             data = st_as_sf(visitDataSpatial),
             options=list(control.fixed = fixedHyper,
                          control.inla= list(int.strategy='eb'),
                          control.compute = list(waic = TRUE, dic = FALSE, cpo = TRUE),
                          verbose = TRUE))

# Assign model summary object and output
modelSummary <- summary(model); modelSummary

# PREDICT -----------------------------------------

# Create grid prediction pixels
ppxl <- mask(UK_R, smoothUK) %>%
  crop(.,smoothUK ) %>%
  as.points %>%
  st_as_sf

# Create multi-time period prediction pixels
ppxlAll <- fm_cprod(ppxl, data.frame( iYear = seq_len(nYear)))

# Predict using spatio-temporal model
# ( N.B. exlcuding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
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
                                                             coverBF +
                                                             coverCF +
                                                             connectivity +
                                                             BFconnINT +
                                                             CFconnINT +
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
linearEffLabels <- c('coverBF' = "Broadleaf cover",
                     'coverCF' = "Conifereous cover",
                     'connectivity' = "Connectivity",
                     'BFconnINT' = "Broadleaf:connectivity",
                     'CFconnINT' = "Conifereous:connectivity",
                     'visitLengthSingle' = "Single record visit",
                     'visitLengthShort' = "Short visit (2-3 records)")
intLabels <- c('BF_pred' = "Broadleaf ",
               'CF_pred' = "Coniferous")

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
  # geom_text(data = timeLabels,
  #           mapping = aes(x = 500, y = 925, label = label)) +
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
  # geom_text(data = timeLabels,
  #           mapping = aes(x = 500, y = 925, label = label)) +
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
  # geom_text(data = timeLabels,
  #           mapping = aes(x = 500, y = 925, label = label)) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        strip.text.x = element_blank()) +
  coord_fixed() +
  geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black", inherit.aes = FALSE)

# COVER-CONNECTIVITY INTERACTION PLOTS

### Set up prediction data frames

# How many prediction steps?
nSamp <- 100

# Extract max scaled value
maxConnectivity <- global(connW, fun = "max", na.rm = TRUE) %>% 
  max

# Create unscaled data frame of cover and connectivity values to predict over
# Separate broadleaf and coniferous data frames
BF_pred_df <-  expand.grid(BF_pred = seq(0, 1, by = 1/nSamp),
                           conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))
CF_pred_df <- expand.grid(CF_pred = seq(0, 1, by = 1/nSamp),
                          conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))

### Scale covariates

# Scale the prediction steps for broadleaf and coniferous woodland separately, and connectivity
# N.B. Have to name columns the same as the original datasets!
BF_pred_df$coverBF_scaled <- ( BF_pred_df$BF_pred - 
                               scalingParams[scalingParams$variable == "coverBF", "variableMean"] ) /
  scalingParams[scalingParams$variable == "coverBF", "variableSD"]

CF_pred_df$coverCF_scaled <- ( CF_pred_df$CF_pred - 
                               scalingParams[scalingParams$variable == "coverCF", "variableMean"] ) /
  scalingParams[scalingParams$variable == "coverCF", "variableSD"]

BF_pred_df$connW_scaled <- CF_pred_df$connW_scaled <- # N.B. Connectivity is the same for both cover types
  (BF_pred_df$conn_pred - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
  scalingParams[scalingParams$variable == "connW", "variableSD"]

# Calculate scaled interaction terms for prediction
BF_pred_df$coverBF_connW <- BF_pred_df$coverBF_scaled * BF_pred_df$connW_scaled
CF_pred_df$coverCF_connW <- CF_pred_df$coverCF_scaled * CF_pred_df$connW_scaled

### Create 'unscaled' vectors of covariate values where species is present for plot

# Subset covarValues to presence records, and then unscale for 
# broadleaf and coniferous woodland, and connectivity
presentCoverBF <- subset(covarValues, presence == 1)$coverBF
presentCoverBF <-
  ((presentCoverBF * scalingParams[scalingParams$variable == "coverBF", "variableSD"]) +
     scalingParams[scalingParams$variable == "coverBF", "variableMean"])

presentCoverCF <- subset(covarValues, presence == 1)$coverCF
presentCoverCF <-
  ((presentCoverCF * scalingParams[scalingParams$variable == "coverCF", "variableSD"]) +
     scalingParams[scalingParams$variable == "coverCF", "variableMean"])

presentConnW <- subset(covarValues, presence == 1)$connW
presentConnW <-
  ((presentConnW * scalingParams[scalingParams$variable == "connW", "variableSD"]) +
     scalingParams[scalingParams$variable == "connW", "variableMean"])

# Join unscaled covariate values together, and take unique values to speed up later steps
presentFixedEff <- data.frame(presentCoverBF, presentCoverCF, presentConnW) %>%
  distinct

# Remove obsolete objects
rm(presentCoverBF, presentCoverCF, presentConnW)

### Predict

# Predict broadleaf cover and connectivity interaction at link scale
BFconnINTpred <- predict(model,
                         BF_pred_df,
                         formula = ~ coverBF_eval(coverBF_scaled) +
                           connectivity_eval(connW_scaled) +
                           BFconnINT_eval(coverBF_connW),
                         exclude = c("spaceTime", "week",
                                     "soilM",  "WMIN", "tasCV", "GDD5", "RAIN",
                                     "coverCF", "CFconnINT"))

# Predict coniferous cover and connectivity interaction at link scale
CFconnINTpred <- predict(model,
                         CF_pred_df,
                         formula = ~ coverCF_eval(coverCF_scaled) +
                           connectivity_eval(connW_scaled) +
                           CFconnINT_eval(coverCF_connW),
                         exclude = c("spaceTime", "week",
                                     "soilM",  "WMIN", "tasCV", "GDD5", "RAIN",
                                     "coverBF", "BFconnINT"))

# Rename prediction columns needed for plot(median)
BFconnINTpred <- rename(BFconnINTpred, BFmedian = median)
CFconnINTpred <- rename(CFconnINTpred, CFmedian = median)

# Join into a single dataframe
allPred <- cbind(BFconnINTpred[, c("BF_pred", "conn_pred", "BFmedian" )], 
                 CFconnINTpred[, c("CF_pred", "CFmedian" )])

### Filter predictions by cover and connectivity values which species is present in

# Find prediction values (discrete) which species fixed effect values overlap with
# Using both broadleaf and connectivity columns of unscaled covariate values where the species is present,
# find if there are any of these values that is within the prediction bin (i.e. x +- max/nSamp)
# for the respective covariate
allPred$BFpresent <- mapply( x = allPred$BF_pred, y = allPred$conn_pred,
                             FUN = function(x, y)
                               
                               any((x - 1 / nSamp) < presentFixedEff$presentCoverBF & 
                                     (x + 1 / nSamp) > presentFixedEff$presentCoverBF  &
                                     (y - maxConnectivity / nSamp) < presentFixedEff$presentConnW  &
                                     (y + maxConnectivity / nSamp) > presentFixedEff$presentConnW ))

allPred$CFpresent <- mapply( x = allPred$CF_pred, y = allPred$conn_pred,
                             FUN = function(x, y)
                               
                               any((x - 1 / nSamp) < presentFixedEff$presentCoverCF &
                                     (x + 1 / nSamp) > presentFixedEff$presentCoverCF  &
                                     (y - maxConnectivity / nSamp) < presentFixedEff$presentConnW  &
                                     (y + maxConnectivity / nSamp) > presentFixedEff$presentConnW ))

# Create separate column where prediction is replaced with 'NA' if no presence records from the cell
allPred <-  allPred %>%
  mutate(BFmedianPres = if_else(BFpresent == TRUE, BFmedian, NA)) %>%
  mutate(CFmedianPres = if_else(CFpresent == TRUE, CFmedian, NA))

# Convert to long data format
allPred <- allPred %>% 
  gather(coverType, cover_pred, "BF_pred" , "CF_pred") %>%
  mutate(median = if_else(coverType == "BF_pred", BFmedian , CFmedian )) %>%
  mutate(medianPres = if_else(coverType == "BF_pred", BFmedianPres , CFmedianPres ))

### Plot

# Plot broadleaf and coniferous cover-connectivity interaction
# (contours based on all predictions, fill only uses prediction space where species is present)
intPlot <- ggplot(allPred, aes(x = cover_pred, y = conn_pred, z = median)) +
  ggtitle("Woodland cover - connectivity interaction") +
  facet_wrap( ~coverType, labeller = as_labeller(intLabels)) +
  geom_tile(aes(fill = medianPres, colour = medianPres)) +
  stat_contour(bins = 50, colour = "black") +
  scale_fill_distiller(na.value = NA,
                       palette = "RdYlBu",
                       direction = 1,
                       guide = guide_colourbar(title = "Posterior\nmedian"),
                       limits = c(-1,1) * max(abs(na.omit(allPred$medianPres)))) +
  scale_colour_distiller(na.value = NA,
                         palette = "RdYlBu",
                         direction = 1,
                         guide = "none",
                         limits = c(-1,1) * max(abs(na.omit(allPred$medianPres)))) +
  scale_x_continuous(name="Proportion woodland cover") +
  scale_y_continuous(name="Connectivity (A)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1),
        panel.grid.major = element_line(colour = "darkgrey"),
        strip.text.x = element_text(size = 12))

### JOINT EVALUATION PLOT

evalPlot <- arrangeGrob(predMedian, predSD ,
                        fixedEffPlot, intPlot ,
                        randomEffPlot, spaceTimePlot,
                        nrow = 3, ncol = 4,
                        widths=c(1, 1, 1, 3), 
                        layout_matrix = rbind(c(1, 1, 1, 2),
                                              c(3, 3, 4, 4),
                                              c(5, 5, 5, 6)),
                        top = grid::textGrob(paste0(iSpecies, ", logCPO = ", logCPO),
                                             gp = grid::gpar(fontsize=20)))

# Compose matern plot
spdeAndMaternPlot <- arrangeGrob(range.plot, covplot,
                                 var.plot, corplot, ncol = 2)

# SAVE -----------------------------------------

# Define species directory
iSpeciesDir <- paste0("Treescapes/ST_SDMs/Output/", 
                      taxaGroup, "/", 
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
       width = 6000, height = 6000, 
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
# # Option 2:
# 
# # Enter directly as component using inbuilt spatial grid data frame recognition
# 
# ... + soilST(main = soilM_SGDF, main_layer = iYear, model = "linear")




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

