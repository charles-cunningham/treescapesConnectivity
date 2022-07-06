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

# LOAD DATA ------------------------------------------

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.prettywkt?download",
              destfile = "Data/Spatial_data/bng.prj")

bng <- "Data/Spatial_data/bng.prj"

# Read in country outline from GADM
UK <- gadm_sf_loadCountries( fileNames = "GBR", level = 1,
                             basefile = "./Data/Spatial_data/")$sf

# Set CRS
UK <- st_transform(UK, st_crs(bng))

# Load in Carabid test data
load("Data/Species_data/Carabids_170316_Cleaned_Data.rdata")

# LOAD FUNCTIONS ------------------------------------

# Load grid coordinate to meters conversion function
source("Functions/gridCoords().R")

# SET PARAMETERS ------------------------------------

# Species to model
mySpecies <- "COL_6133"

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 80 ### Set v coarse for testing, best value will probably be ~50

### CODE --------------------------------------------

### PROCESS GB SPATIAL FILES ------------------------

# Select countries

# Remove Northern Ireland (not enough data?)
# GB <- test[test$NAME_1 != "Northern Ireland", "NAME_1"]

# Only Wales for testing
GB <- test[test$NAME_1 == "Wales", "NAME_1"]

# Write layer to shapefile
st_write(GB, "Data/Spatial_data/GB.shp", append = FALSE)

# Read in as SpatVector
GB <- vect("Data/Spatial_data/GB.shp")

# Disaggregate
GB <- disagg(GB)

# Calculate area
GB$area_sqkm <- expanse(GB, unit="km")

# Remove polygons with < 20km ^2 area
GB <- GB[GB$area_sqkm > 10]

# Aggreagte back
GB <- aggregate(GB) #aggregate back

# # Create base grid (units are metres)
# baseR <- rast(crs = bng, 
#              xmin = 0, xmax = 700000,
#              ymin = 0, ymax = 1300000,
#              resolution = 1000)
# 
# # Rasterise, including cells that touch the edge
# GBr <- rasterize(x = GB, y = baseR, touches = TRUE)
# 
# plot(GBr)

### PROCESS TO VISITS DATAFRAME -----------------------

# Add visit column (gridref plus date)
taxa_data$visit <- paste(taxa_data$TO_GRIDREF, taxa_data$TO_STARTDATE, sep="_")

# Transform into visit data.frame
visitData <- taxa_data %>%
  add_column(Presence = 1) %>% # Add a presence column
  pivot_wider(names_from = CONCEPT,
              values_from = Presence) %>% # Transform to wide data frame (species names now columns)
  janitor::clean_names(case = "all_caps") %>% # All characters in data.frame to CAPS
  rename(Presence = toupper(mySpecies)) %>% # Change the species of interest column name to 'Presence'
  dplyr::select(TO_GRIDREF,TO_STARTDATE,
                YEAR, VISIT, Presence) # Convert back to long data.frame

# Change non-detections visits to zeros
visitData$Presence[is.na(visitData$Presence)] <- 0

# Check number of presense and pseudo-absences
table(visitData$Presence)

# CREATE EFFORT COVARIATES --------------------------

### Add day of year column ( will be included in model as covariate)
visitData$WOY <- visitData$TO_STARTDATE %>%
  strftime(., format = "%V") %>%
  as.numeric(.)

### Calculate "flight period", i.e. approximate dates when species is observable

# Create empty vector to hold flight period days of year
flightPeriodAllWeeks <- c()

# Loop through visits where species was detected...
for ( i in visitData$WOY[ visitData$Presence == 1 ]) {
  
  # Save all days +/- 5 from positive records
  flightPeriodAllWeeks <- c(flightPeriodAllWeeks, (i-1) : (i+1))

}

# Correct weeks of year that have fallen into other years
 flightPeriodAllWeeks[flightPeriodAllWeeks < 1] <- flightPeriodAllWeeks[flightPeriodAllWeeks < 1] + 52
 flightPeriodAllWeeks[flightPeriodAllWeeks > 52] <- flightPeriodAllWeeks[flightPeriodAllWeeks > 52] - 52

# Get unique values
flightPeriodAllWeeks <- unique(flightPeriodAllWeeks) 

### For each unique location, find number of visits in flight period
# (May be faster way to do this)

# Add number of visits column
visitData$nVisits <- 0

# Loop through unique GRIDREF values (locations)
for (i in unique(visitData$TO_GRIDREF)) {
  
  # Find number of visits in flight period
  nVisits <- visitData[visitData$TO_GRIDREF == i & #Subset data to GRIDREF location, plus ...
                          visitData$WOY %in% flightPeriodAllWeeks,] %>% # ... day of year in flight period days, 
    NROW(.) #  and find number of rows (visits)

  # Assign nVisits to column of location
  visitData[visitData$TO_GRIDREF == i,]$nVisits <- nVisits
  
  }

# - Effort: Add log of number of visits to site in year as measure of effort 
# (but only in flight period, i.e. within +/- 5 days of any positive record)?

# CONVERT GRIDREF TO METERS -------------------------

# Create x and y columns
visitData$X <- NA
visitData$Y <- NA

# Populate x and y columns using gridCoords function
# N.B. Function find the bottom, left-hand corner of the 1km grid square, 
# so need to add 500m

for (i in 1:NROW(visitData)) { # Loop through each row
  
  # Add in x and y
  visitData[i, "X"] <- gridCoords(visitData[i, "TO_GRIDREF"], units = "m")$x + 500
  visitData[i, "Y"] <- gridCoords(visitData[i, "TO_GRIDREF"], units = "m")$y + 500
  
}

### PROCESS COVARIATES -----------------------

# Scale covariates
visitData$WOYscaled <- as.numeric(scale(visitData$WOY))
visitData$nVisitsLog <- as.numeric(log(visitData$nVisits))

### CONVERT TO SPATIAL DATA

# Make sf object
visitDataSpatial <- visitData %>%
  st_as_sf(.,coords=c("X","Y"), crs = bng)

# Plot
ggplot() +
  geom_point(data = subset(visitData, Presence==0), aes(x=X, y=Y),
             color = "white") +
  geom_point(data = subset(visitData, Presence==1), aes(x=X, y=Y),
             color = "black") +
  geom_sf(data = st_as_sf(GB), colour = "black", fill = NA)

# Some coordinates in the Irish sea so need to filter these

# Filter to GB
visitData <- visitData %>%
  filter(st_intersects(visitDataSpatial, st_as_sf(GB), sparse = FALSE)[,1]) 


### CONVERT TO 'SP' AND KM RESOLUTION FOR INLABRU (CUE WARNING MESSAGES)

visitData_sp <- visitData
coordinates(visitData_sp) <- c("X", "Y")
proj4string(visitData_sp) <- st_crs(bng)$proj4string
GB_sp <- as(st_as_sf(GB), 'Spatial')
proj4string(GB_sp) <- st_crs(bng)$proj4string

# Reproject to km for inlabru (needs changing back after!)
GB_sp <- spTransform(GB_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))
visitData_sp <- spTransform(visitData_sp, gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# Only  visitData_sp and GB_sp converted for inlabru analysis below

### CREATE MESH --------------------------------

# Subset to several years for testing
visitData_sp <- visitData_sp[ visitData_sp$YEAR > 1985 & visitData_sp$YEAR <= 1990 ,]

# Max edge is range/8 as a rule of thumb
maxEdge = estimated_range/8

mesh <- inla.mesh.2d(boundary = GB_sp,
                     loc = visitData_sp,
                     max.edge = c(1,4) * maxEdge,
                     cutoff = maxEdge/2,
                     offset =c(1,4) * maxEdge,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

ggplot() + 
  gg(mesh) + 
  geom_sf(data = st_as_sf(subset(visitData_sp, Presence==0)), color = "white") +
  geom_sf(data = st_as_sf(subset(visitData_sp, Presence==1)), color = "black") +
  coord_fixed() +
  geom_sf(data = st_as_sf(GB_sp), col = "black", fill = NA)

### FIT SPATIO-TEMPORAL MODEL ----------------

# Simplify year to nearest decade (trying annual model now so redundant)

# floor_decade <- function(value){ return(value - value %% 10) }
# visitData_sp$decade <- floor_decade((visitData_sp$YEAR))
# visitData_sp$iYear <- ((visitData_sp$decade - min(visitData_sp$decade)) / 10) + 1

# Create indices
visitData_sp$iYear <- visitData_sp$YEAR - min(visitData_sp$YEAR) + 1
iYear <- visitData_sp$iYear
nYear <- length(unique(visitData_sp$iYear))

# Define spatial pattern
mySpace <- inla.spde2.pcmatern(
  mesh,
  alpha = 2,
  prior.range = c(10 * maxEdge, 0.5),   
  prior.sigma = c(1, 0.5))

### Simple space-time model with effort (linear log(nVisits)) and spatial field

# Specify components
cmp <- Presence ~  nVisitsLog +
  spatial(main = coordinates,
          group = iYear,
          ngroup = nYear,
          model = mySpace,
          control.group=list(model="ar1")) +
  Intercept(1)

# Fit model
inlabruST <- bru(cmp,
                family = "binomial",
                data = visitData_sp,
                options=list(verbose=TRUE))

# Model summary
summary(inlabruST)

### Trialling week of year with smoothing process, currently rw2 but perhaps splines would be best?

inlabruCmp  <-  ~ 
  WOY(main = WOY , model = "rw2") +
  nVisitsLog(main = nVisitsLog , model = "linear") + 
  spatial(main = coordinates,
          group = iYear,
          ngroup = nYear,
          model = mySpace,
          control.group=list(model="ar1")) +
  Intercept(1)

inlabruFormula <- Presence ~ 
  spatial + Intercept + WOY + nVisitsLog
  
# Fit model
inlabruST2 <- bru(components = inlabruCmp,
                 formula = inlabruFormula,
                 family = "binomial",
                 data = visitData_sp,
                 options=list(verbose=TRUE))

summary(inlabruST2)

### PREDICT -----------------------------------------

# Create spatial data frame to predict over
ppxl <- pixels(mesh, mask = GB_sp)
ppxlAll <- cprod(ppxl, data.frame(iYear = seq_len(nYear)))

# Predict using spatio-temporal model
inlabruPred <- predict(inlabruST2, 
                       ppxlAll, 
                       ~ data.frame(iYear = iYear,
                                    lambda =  plogis(Intercept + 
                                                       spatial +
                                                       WOY +
                                                       nVisitsLog)))

plot(inlabruST2, "WOY")

# Plot
ggplot() + 
  gg(inlabruPred, aes(fill=mean)) +
  scale_fill_viridis_c(option="magma",direction=1) +
  geom_sf(data = st_as_sf(subset(visitData_sp, Presence==1)), color = "blue") +
  geom_sf(data = st_as_sf(subset(visitData_sp, Presence==0)), color = "red") +
  facet_wrap(~iYear) +
  coord_fixed() +
  geom_sf(data = st_as_sf(GB_sp), col = "black", fill = NA)
