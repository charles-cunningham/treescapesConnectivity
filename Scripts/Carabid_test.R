# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: 
#
# Script Description: 
#
#
# Notes: Trying to work in 'sf' as much as possible and remove 'sp'
# i.e. no 'raster' package!
#

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(INLA) 
library(inlabru)
library(GADMTools)
library(terra)
library(tidyverse)

# LOAD DATA AND FUNCTIONS ---------------------------

# Download BNG WKT string
download.file(url = "https://epsg.io/27700.prettywkt?download",
              destfile = "Data/Spatial_data/bng.prj")

bng <- "Data/Spatial_data/bng.prj"

# Load grid coordinate to meters conversion function
source("Functions/gridCoords().R")

# Load in Carabid test data
load("Data/Species_data/Carabids_170316_Cleaned_Data.rdata")

# SET PARAMETERS ------------------------------------


### CODE --------------------------------------------

### CREATE GB SPATIAL FILE --------------------------

# Read in country outline from GADM
gadmUK <- gadm_sf_loadCountries( fileNames = "GBR", level = 1, basefile = "./Data/Spatial_data/")

# Remove N. Ireland as no data for here unfortunately
gadmGB <- gadm_remove(gadmUK, regions="Northern Ireland", 1)

# Save to shapefile, if file has not been created yet
if (!file.exists("Data/Spatial_data/GB.shp")) {
gadm_exportToShapefile(gadmGB, "Data/Spatial_data/GB.shp") }

# Read in as SpatVector
GB <- vect("Data/Spatial_data/GB.shp")

# Reproject to BNG
GB <- project(GB, bng)

# Disaggregate
GB <- disagg(GB)

# Calculate area
GB$area_sqkm <- expanse(GB, unit="km")

# Remove polygons with < 20km ^2 area
GB <- GB[GB$area_sqkm > 10]

# Aggreagte back
GB <- aggregate(GB) #aggregate back

# Create base grid (units are metres)
baseR <- rast(crs = bng, 
             xmin = 0, xmax = 700000,
             ymin = 0, ymax = 1300000,
             resolution = 1000)

# Rasterise, including cells that just touch the edge
GBr <- rasterize(x = GB, y = baseR, touches = TRUE)

#plot(GBr)

### PROCESS TO USABLE DATAFRAME ---------------------

# CONVERT GRIDREF TO METERS -------------------------

# Create x and y columns
taxa_data$X <- NA
taxa_data$Y <- NA

# Populate x and y columns using gridCoords function
# N.B. Function find the bottom, left-hand corner of the 1km grid square, 
# so need to add 500m

for (i in 1:NROW(taxa_data)) { # Loop through each row
  
  # Add in x and y
  taxa_data[i, "X"] <- gridCoords(taxa_data[i, "TO_GRIDREF"], units = "m")$x + 500
  taxa_data[i, "Y"] <- gridCoords(taxa_data[i, "TO_GRIDREF"], units = "m")$y + 500
  
}

# KEEP USEFUL COLUMNS
taxa_data <- taxa_data[,c("CONCEPT","YEAR","X", "Y")] 


# Get numbers of unique values
nRecords <- data.frame("CONCEPT" = unique(taxa_data["CONCEPT"]),
                       "n" = NA)

for (i in 1:NROW(nRecords)) {
  
  nRecords[i,"n"] <- taxa_data[nRecords[i,"CONCEPT"] == taxa_data$CONCEPT,] %>% 
    NROW(.)
   
  
}

nRecords[nRecords$n ==  max(nRecords$n),]

# Col_6133

testData <- taxa_data[taxa_data$CONCEPT == "Col_6133" ,]

test <- vect(testData, geom=c("X", "Y"), crs= bng) %>%
  rasterize(., GBr)
plot(test, add = TRUE, col = "black")


# have test and GBr, create data frame for use in unmarked

### UNMARKED TEST

# Create data frame of all coords
testDF <- crds(GBr, df=TRUE, na.rm=TRUE)

# Add each year as a column
testDF[c(min(testData$YEAR): max(testData$YEAR)) %>%
  as.character(.) ] <- 0

# Assign records from testdata to testDF
for (i in 1:NROW(testData)) {
  
  # Extract year and coords from from of testdata
  year <- testData$YEAR[i] %>% as.character(.)
  coords <- c(testData$X[i], testData$Y[i])
  
  # Change value of testDF to 1 to indicate record
  testDF[testDF$x == coords[1] & testDF$y == coords[2],
         year] <- 1

}
siteCovs <-  testDF[ ,1:2]


umf <- unmarkedFrameOccu( y = testDF[,3:NCOL(testDF)] , siteCovs = siteCovs)

summary(umf)

fm <- occu(formula = ~ 1 
           ~ x + y,
           data = umf)

fm

backTransform(fm, type = "state")

occuPred <- predict(fm,
                    type = "state",
                    newdata = testDF[ ,1:2],
                    na.rm = TRUE,
                    inf.rm = TRUE)


levelplot(Predicted ~ testDF$x + testDF$y,
          data = occuPred,
          col.regions = rev(terrain.colors(100)),
          at = seq(0,1,length.out=101))


### SPOCCUPANCY TEST

data(hbef2015)
sp.names <- dimnames(hbef2015$y)[[1]]
btbwHBEF <- hbef2015
btbwHBEF$y <- btbwHBEF$y[sp.names == "BTBW", , ]


str(hbef2015)

(hbef2015$y)

matrix()

# Create data frame of all coords
coords <- crds(GBr, df=FALSE, na.rm=TRUE)

str(testDF)

# Add each year as a column
testDF[c(min(testData$YEAR): max(testData$YEAR)) %>%
         as.character(.) ] <- 0

# Assign records from testdata to testDF
for (i in 1:NROW(testData)) {
  
  # Extract year and coords from from of testdata
  year <- testData$YEAR[i] %>% as.character(.)
  coords <- c(testData$X[i], testData$Y[i])
  
  # Change value of testDF to 1 to indicate record
  testDF[testDF$x == coords[1] & testDF$y == coords[2],
         year] <- 1
  
}




