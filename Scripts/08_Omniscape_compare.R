# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Analyse Omniscape output 
#
# Script Description: Compare the effect of different parameters
# on the connectivity of different landscapes

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(tidyverse)
library(terra)
library(corrplot)

# Read in LCM SpatRasters
UKLCM1990 <- rast("../Data/Spatial_data/LCM/1990/UK1990all_1km.tif")[["W"]]
UKLCM2015 <- rast("../Data/Spatial_data/LCM/2015/UK2015all_1km.tif")[["W"]]

# Change extent (otherwise western UK-Ireland border too close to edge of extent)
UKLCM1990 <- extend(UKLCM1990, ext( -10000, 700000, 0, 1300000))
UKLCM2015 <- extend(UKLCM2015, ext( -10000, 700000, 0, 1300000))

# Change names
names(UKLCM1990) <- "Cover1990"
names(UKLCM2015) <- "Cover2015"

# CREATE CONNECTIVITY TIFS --------------------------

# Read input table
inputTable <- read.table("../Data/Spatial_data/Omniscape/inputTable.txt",
                         header = TRUE)

omniOutput <- lapply (1:NROW(inputTable), function(x) { ####!!! Change this line when final connectivity run finishes
  
  # Create raster fro row i of input table
  omniR <- paste0(
    "../Data/Spatial_data/Omniscape/",
    "radius",
    inputTable$Radius[x],
    "_resistance",
    inputTable$Resistance[x],
    "/output",
    inputTable$Year[x],
    "/cum_currmap.tif"
  ) %>%
    rast(.)
  
  # Set rast names
  names(omniR) <- paste0("radius", inputTable$Radius[x],
                         "_resistance", inputTable$Resistance[x],
                         "_year", inputTable$Year[x])
  return(omniR)
  
}) %>% rast(.)

# AGGREGATE AND RESCALE

# Aggregate to 1km
omniOutput1km <- terra::aggregate(omniOutput, fact = 40)

# Rescale
#omniOutput1km <- terra::scale(omniOutput1km)

# COMPARE LCM WITH OMNISCAPE OUTPUT --------------------------

# Create single data frame with LCm and omniscape output
all_df <- c( UKLCM1990,  UKLCM2015, omniOutput1km) %>%
  as.data.frame() %>%
  drop_na

# Create data frames for each time period
y1990_df <- all_df[,grep("1990", c( names(UKLCM1990), names(UKLCM2015), names(omniOutput1km)) )] 
y2015_df <- all_df[,grep("2015", c( names(UKLCM1990), names(UKLCM2015), names(omniOutput1km))  )]

# CORRELATION

res <- cor(y1990_df)
#png(filename = 'figures/corr_plot.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res, type = "upper", order = "hclust",
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)

res <- cor(y2015_df)
#png(filename = 'figures/corr_plot.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res, type = "upper", order = "hclust",
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)


### testing
# 
# test <- data.frame(UKLCM1990[], omniOutput1km[[29]][]) %>% na.omit(.)
# 
# ggplot(test, aes(x = W, y = test[,2])) +
#   geom_point(size = 0.1)
#   
# 
# 
# plot(UKLCM2015[], omniOutput1km[[10]][],
#      xlab = "Proportion 1km woodland cover",
#      ylab = "Normalised 1km connectivity (mean of pixels). Radius: 40. Resistance: 100")
# 
# 
# plot(omniOutput1km[[10]])
# 
# y1990_df[,1] %>%
#   log(.) %>%
#   .[is.finite(.)] %>%
#   scale(.) %>%
#   hist(.)
# 
# omniOutput1km[[7]][] %>% hist(.)
# 
# test <- omniOutput1km[[7]]
# 
