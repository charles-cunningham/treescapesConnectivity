# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Process LCM and create files for Omniscape run  
#
# Script Description: This script (1) processes LCM woodland 
# cover maps for use by Omniscape, and (2) creates directories,
# .INI files etc that are required for Omniscape to run.

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)

# Set fraction of RAM that may be used by the program
terraOptions(memfrac = 0.9)

# PARAMETERS ------------------------------------------

# Specify parameter options for radius and resistance here
radius     <- c(20, 40, 80, 160) # 25m pixels, so 20 = 500m, 40 = 1km, 80 = 2km, 160 = 4km
blockSize  <- c(1, 3, 5, 15) # Rule of Thumb : must be <10% of radius
year       <- c(1990, 2015)
resistance <- c(10, 100, 1000)

# Set arbitrary reassignment value
arbitraryValue <- 2

# DATA FILES ------------------------------------------

### Download BNG WKT string
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = "Treescapes/Omniscape/bng.prj")

bng <- sf::st_crs("Treescapes/Omniscape/bng.prj")$wkt

### Land Cover Maps of woodland

# 1990
LCM1990 <- rast("Treescapes/Omniscape/UK1990all.tif")
# 2015
LCM2015 <- rast("Treescapes/Omniscape/UK2015all.tif")

# # 1990
# LCM1990 <- rast("../Data/Spatial_data/LCM/1990/UK1990all.tif")
# # 2015
# LCM2015 <- rast("../Data/Spatial_data/LCM/2015/UK2015all.tif")

### National boundaries

# UK
UK <- vect("Treescapes/Omniscape/GBR.shp")
# Ireland
Ireland <- vect("Treescapes/Omniscape/IRL.shp")

# # UK 
# UK <- vect("../Data/Spatial_data/Boundaries_and_CRS/UK/GBR.shp")
# # Ireland 
# Ireland <- vect("../Data/Spatial_data/Boundaries_and_CRS/Ireland/IRL.shp")

### PROCESS LAND COVER MAPS FOR OMNISCAPE ----------------------

# Select all woodland layer
omni1990W <- LCM1990[["W"]]
omni2015W <- LCM2015[["W"]]

### Add in land from Ireland

# Change extent (otherwise western UK-Ireland border too close to edge of extent)
omni1990W <- extend(omni1990W, ext( -10000, 700000, 0, 1300000))
omni2015W <- extend(omni2015W, ext( -10000, 700000, 0, 1300000))

# Create 'background' raster *including Ireland* with 0 values
allLand <- rbind(UK, Ireland) %>% # Join together spatVectors
  aggregate(.) %>% # Remove boundaries
  rasterize(., omni1990W, 0, touches  = TRUE) # Rasterise using woodland spatVector

# Join background and UK LCM map (i.e. include Ireland)
omni1990W <- mosaic(omni1990W, allLand, fun = "max")
omni2015W <- mosaic(omni2015W, allLand, fun = "max")

# Reclass '0' (non-woodland terrestrial cells) to an arbitrary >0 integer
# (Omniscape cannot accept 0 values, also can't assign '1' as that's woodland!)
omni1990W <- subst(omni1990W, from = 0, to = arbitraryValue)
omni2015W <- subst(omni2015W, from = 0, to = arbitraryValue)

# Save for OMNISCAPE run
writeRaster(omni1990W,
            "Treescapes/Omniscape/omni1990Cover.tif",
            overwrite = TRUE)
writeRaster(omni2015W,
            "Treescapes/Omniscape/omni2015Cover.tif",
            overwrite = TRUE)

# CREATE PARAMETERISATION TABLE ------------------------------

# Create input table to populate different parameter iterations
inputTable <- data.frame( "Radius" = rep(radius,
                                         each = length(year) *
                                           length(resistance)),
                          "BlockSize" = rep(blockSize,
                                            each = length(year) *
                                              length(resistance)),
                          "Resistance" = rep(resistance,
                                             each = length(year)),
                          "Year" = year)

# Write input table to tab-seperated .txt file
write.table(inputTable,
            paste0("Treescapes/Omniscape/inputTable.txt"),
            sep = "\t",
            row.names = FALSE)

# CREATE OMNISCAPE FILES ------------------------------------

# Loop though every row of input table to populate Omniscape folders
for (i in 1:NROW(inputTable)) {

  ### Create directory

  # Specify omniscape folder string for row 'i'
  iDir <-
    paste0(
      "Treescapes/Omniscape/radius",
      inputTable$Radius[i],
      "_resistance",
      inputTable$Resistance[i]
    )

  # Create iDir
  dir.create(iDir,
             showWarnings = FALSE)

  ### Create omniscape reclassification table within folder 'i'

  # Create reclass table using resistance 'i'
  reclassTable <-
    matrix(c(arbitraryValue,
             inputTable$Resistance[i]),
           ncol = 2)

  # Write reclass table to tab-seperated .txt file
  write.table(reclassTable,
              paste0(iDir, "/reclassTable.txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)

  ### Create omniscape .INI file for iteration 'i' (resistance not set here)
  # N.B. Full file path needed

  textINI <- paste0( "[Required]

; Set treescape layer to use #!# Changes #!#
resistance_file = ", "/users/cac567/scratch/Treescapes/Omniscape/omni", inputTable$Year[i], "Cover.tif", "

; Set radius in pixels (25m) #!# Changes #!#
radius = ", inputTable$Radius[i], "

; Set output folder name #!# Changes #!#
project_name = ", "/users/cac567/scratch/", iDir, "/output", inputTable$Year[i], "

[General options]

; Block size #!# Changes #!#
block_size = ", inputTable$BlockSize[i], "

; Should a source layer be derived using the resistance layer?
source_from_resistance = true

; Maximum resistance value a cell can have to be included as a source
r_cutoff = 1

; Solver (cholmod seems to perform better)
solver = cholmod

[Resistance Reclassification]

; Reclassify resistance?
reclassify_resistance = true

; File path to reclass table #!# Changes #!#
reclass_table = ", "/users/cac567/scratch/", iDir, "/reclassTable.txt", "

; Save reclassified raster to check?
write_reclassified_resistance = true

[Processing options]

; Parallelize?
parallelize = true

; Batch size
parallel_batch_size = 10

[Output options]

; Write raw cumulative current to disk?
write_raw_currmap = true

; Specify whether to mask current flow output with resistance NoData values
mask_nodata = true")

  # Save to .INI file
  writeLines(textINI,
             paste0( iDir,
                     "/omniscapeSettings", inputTable$Year[i], ".ini"))  #!# Changes #!#
}

# CALCULATE JOB SUBMISSION ARRAY SIZES ------------------------

# Total array size
print( paste("The total array size is",
             NROW(inputTable) ))

# Array sizes for radii
# N.B. Need to submit separately as vastly different computing times

for (i in unique(inputTable$Radius)) {

  print( paste0("The array for radius ",
               i,
               " is ",
               min(which(inputTable$Radius == i)),
               "-",
               max(which(inputTable$Radius == i))))
}
