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

# # CREATE CONNECTIVITY TIFS --------------------------
# # Takes ~30 mins to run so save at end of section to only need to run once
#
# # Read input table
# inputTable <- read.table("../Data/Spatial_data/Omniscape/inputTable.txt",
#                          header = TRUE)
# 
# omniOutput <- lapply (1:NROW(inputTable), function(x) { ####!!! Change this line when final connectivity run finishes
#   
#   # Create raster for row i of input table
#   omniR <- paste0(
#     "../Data/Spatial_data/Omniscape/",
#     "radius",
#     inputTable$Radius[x],
#     "_resistance",
#     inputTable$Resistance[x],
#     "/output",
#     inputTable$Year[x],
#     "/cum_currmap.tif"
#   ) %>%
#     rast(.)
#   
#   # Set rast names
#   names(omniR) <- paste0("radius", inputTable$Radius[x],
#                          "_resistance", inputTable$Resistance[x],
#                          "_year", inputTable$Year[x])
#   return(omniR)
#   
# }) %>% rast(.)
# 
# # Aggregate to 1km
# omniOutput1km <- terra::aggregate(omniOutput, fact = 40)
# 
# # Rescale
# #omniOutput1km <- terra::scale(omniOutput1km)
# 
# # Save
# writeRaster(omniOutput1km,
#             "../Data/Spatial_data/Omniscape/omniConn_1km.tif",
#             overwrite = TRUE) 

# COMPARE LCM WITH OMNISCAPE OUTPUT --------------------------

# load data
omniOutput1km <- rast("../Data/Spatial_data/Omniscape/omniConn_1km.tif")

# Create single data frame with LCM cover and Omniscape connectivityoutput
all_df <- c( UKLCM1990,  UKLCM2015, omniOutput1km) %>%
  as.data.frame() %>%
  drop_na

# Create data frames for each time period
y1990_df <- all_df[,grep("1990", c( names(UKLCM1990), names(UKLCM2015), names(omniOutput1km)) )] 
y2015_df <- all_df[,grep("2015", c( names(UKLCM1990), names(UKLCM2015), names(omniOutput1km))  )]

# OMNISCAPE COMPARISON PLOT ----------------------------------

# Edit column names for plot
names(y1990_df) <- names(y1990_df) %>%
  gsub("radius", "radius=", .) %>%
  gsub("resistance", "resistance=", .) %>%
  gsub("_", ", ", .) %>%
  gsub(", year1990", "", .) %>%
  gsub("Cover1990", "Cover", .)
names(y2015_df) <- names(y2015_df) %>%
  gsub("radius", "radius=", .) %>%
  gsub("resistance", "resistance=", .) %>%
  gsub("_", ", ", .) %>%
  gsub(", year2015", "", .) %>%
  gsub("Cover2015", "Cover", .)

# CORRELATION PLOT -------------------------------------------------

# 1990 correlation plot
cor1990 <- cor(y1990_df)
png(filename = '../Writing/Plots/corr_omni_1990.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor1990, type = "upper", order = "original",
         method = "square", addCoef.col="white", tl.col = "black",
          tl.srt = 45)
dev.off()

# 2015 correlation plot
cor2015 <- cor(y2015_df)
png(filename = '../Writing/Plots/corr_omni_2015.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor2015, type = "upper", order = "original",
         method = "square", addCoef.col="white",  tl.col = "black",
          tl.srt = 45)
dev.off()

# Combine into single plot
corBoth <- cowplot::ggdraw(clip = "on") +
  cowplot::draw_image('../Writing/Plots/corr_omni_1990.png', -0.05, 0.45, 1.1, 0.55) +
  cowplot::draw_image('../Writing/Plots/corr_omni_2015.png', -0.05, -0.05, 1.1, 0.55) +
  cowplot::draw_label("(a) 1990", 0.1, 0.95, size = 22) +
  cowplot::draw_label("(b) 2015", 0.1, 0.45, size = 22)
  
# Save
ggsave(filename = paste0("../Writing/Plots/", "omniCorr.png"),
       corBoth,
           dpi = 600,
           units = "px", width = 4000, height = 8000)  

# Remove individual plots as not needed
unlink('../Writing/Plots/corr_omni_1990.png')
unlink('../Writing/Plots/corr_omni_2015.png') 

# N.B.
# Radius 160 (4km) and resistance 100 is compromise between correlation
# and previous estimates of conductivity from literature: we use this.
