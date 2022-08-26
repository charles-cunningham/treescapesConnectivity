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
devtools::install_github('BiologicalRecordsCentre/sparta')
library(sparta)
# LOAD DATA ------------------------------------------

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

climateDir <- "../Data/Spatial_data/CHESS/Current_climate_variables"

# LOAD FUNCTIONS ------------------------------------

# Load grid coordinate to meters conversion function
source("Functions/gridCoords().R")

# SET PARAMETERS ------------------------------------

# Species to model
mySpecies <- "COL_6133"

# Estimated range of spatial effect in km (determines mesh)
estimated_range <- 50 ### Set v coarse for testing, best value will probably be ~50

### CODE --------------------------------------------

### PROCESS GB SPATIAL FILES ------------------------

# Select countries

# Remove Northern Ireland (not enough data?)
# GB <- test[test$NAME_1 != "Northern Ireland", "NAME_1"]

# Only Wales for testing
GB <- UK[UK$NAME_1 == "Wales", "NAME_1"]

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

# Aggregate back
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

# CONVERT GRIDREF TO METERS -------------------------

# Create x and y columns
taxa_data$X <- NA
taxa_data$Y <- NA

# Populate x and y columns using gridCoords function
# N.B. Function find the bottom, left-hand corner of the 1km grid square, 
# so need to add 500m

system.time(
for (i in unique(taxa_data$TO_GRIDREF)) { # Loop through each row
   i = unique(taxa_data$TO_GRIDREF)[1]
  # Run gridCoords()
  coordData <- gridCoords(i, units = "m")
  
  # Add in x and y
  taxa_data[taxa_data$TO_GRIDREF == i, "X"] <- coordData$x + 500
  taxa_data[taxa_data$TO_GRIDREF == i, "Y"] <- coordData$y + 500
  
}
)

system.time(
  
  sparta::gr2gps_latlon(taxa_data$TO_GRIDREF[1], precision = 1000, projection = "OSGB", centre = TRUE)
  
)


### PROCESS TO VISITS DATAFRAME -----------------------

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

# CREATE EFFORT COVARIATES --------------------------

### Add day of year column ( will be included in model as covariate)
visitData$WOY <- visitData$TO_STARTDATE %>%
  strftime(., format = "%V") %>%
  as.numeric(.)

# CREATE CLIMATE COVARIATES --------------------------

climateDir



### CONVERT TO SPATIAL DATA

# Make sf object
visitDataSpatial <- visitData %>%
  st_as_sf(.,coords=c("X","Y"), crs = bng)

# Plot
ggplot() +
  geom_point(data = subset(visitData, PRESENCE == 0), aes(x=X, y=Y),
             color = "white") +
  geom_point(data = subset(visitData, PRESENCE == 1), aes(x=X, y=Y),
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
#visitData_sp <- visitData_sp[ visitData_sp$YEAR > 1985 & visitData_sp$YEAR <= 1989 ,]

# Max edge is range/8 as a rule of thumb
maxEdge = estimated_range/8

mesh <- inla.mesh.2d(boundary = GB_sp,
                     loc = visitData_sp,
                     max.edge = c(1,4) * maxEdge,
                     cutoff = maxEdge/2,
                     offset = c(1,4) * maxEdge,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

ggplot() + 
  gg(mesh) + 
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE == 0)), color = "white") +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE == 1)), color = "black") +
  coord_fixed() +
  geom_sf(data = st_as_sf(GB_sp), col = "black", fill = NA)

### PROCESS COVARIATES -----------------------

# Scale covariates
#visitData$WOYscaled <- as.numeric(scale(visitData$WOY))
#visitData$nVisitsLog <- as.numeric(log(visitData$nVisits))

# Sort factor levels (need to convert to dummy variables)
visitData_sp <- visitData_sp %>%
  model.matrix(object = ~VISIT_LENGTH) %>%
  as.data.frame() %>%
  select(-1) %>%
  cbind(visitData_sp, .)


### FIT SPATIO-TEMPORAL MODEL ----------------

# Simplify year to nearest decade (trying annual model now so redundant)

floor_decade <- function(value){ return(value - value %% 10) }
visitData_sp$decade <- floor_decade((visitData_sp$YEAR))
visitData_sp$iYear <- ((visitData_sp$decade - min(visitData_sp$decade)) / 10) + 1

# Create indices
#visitData_sp$iYear <- visitData_sp$YEAR - min(visitData_sp$YEAR) + 1
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
  WOY(main = WOY , model = "rw2") + 
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

summary(inlabruST2)


### PREDICT -----------------------------------------

# Create spatial data frame to predict over
ppxl <- pixels(mesh, mask = GB_sp)
ppxlAll <- cprod(ppxl, data.frame(iYear = seq_len(nYear)))

# Predict using spatio-temporal model
# ( N.B. exlcuding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
inlabruPred <- predict(inlabruST2, 
                       ppxlAll, 
                       ~ data.frame(iYear = iYear,
                                    lambda =  plogis(Intercept + 
                                                       spatial +
                                                       max(inlabruST2$summary.random$WOY$mean)))) # max value for WOY to predict over
# Control for WOY and list length
# Plot
ggplot() + 
  gg(inlabruPred, aes(fill=mean)) +
  scale_fill_viridis_c(option="magma",direction=1) +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE==1)), color = "blue") +
  geom_sf(data = st_as_sf(subset(visitData_sp, PRESENCE==0)), color = "red") +
  facet_wrap(~iYear) +
  coord_fixed() +
  geom_sf(data = st_as_sf(GB_sp), col = "black", fill = NA)
