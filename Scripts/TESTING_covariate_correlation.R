# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Check covariate correlation
#
# Script Description: 

# LOAD LIBRARIES & INSTALL PACKAGES -----------------

# Change  library to C: (R: doesn't have enough space for packages):
.libPaths("R:/rsrch/cb751/lab/Charles/R/PackageLibrary")

# Load packages
library(terra)
library(tidyverse)
library(corrplot)

# DATA FILES ------------------------------------------

# SET PARAMETERS ------------------------------------

# List all covariate spatial grid data frames
# SpatRasters
for (i in list.files("../Data/Spatial_data/DataForInlabru/spatRaster",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("../Data/Spatial_data/DataForInlabru/spatRaster/",
                     i)))
}

cov_R <- list(coverBF_scaled, coverCF_scaled, connW_scaled,
              coverBF_connW, coverCF_connW,
              GDD5_grp, WMIN_grp, tasCV_grp, soilM_grp, RAIN_grp)

covNames <- c("Broadleaf Cover", "Coniferous Cover", "Connectivity",
               "Broadleaf-Connectivity interaction",
               "Coniferous-Connectivity interaction",
               "GDD5", "WMIN", "tasCV", "soilM", "RAIN")

# PROCESS DATA --------------------------------------

cov_R_1990 <- lapply(cov_R, function(x) { x[[1]] }) %>% rast(.)
cov_R_2015 <- lapply(cov_R, function(x) { x[[2]] }) %>% rast(.)

names(cov_R_1990) <- names(cov_R_2015) <- covNames

# CHECK CORRELATION ---------------------------------

full_df <- as.data.frame(cov_R_2015) %>% drop_na(.)
res <- cor(full_df)
#png(filename = 'figures/corr_plot.png', width = 45, height = 30, units = "cm", res = 300)
corrplot(res, type = "upper", order = "hclust",
         method = "color", addCoef.col="black", number.cex=0.75,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8)
dev.off()
