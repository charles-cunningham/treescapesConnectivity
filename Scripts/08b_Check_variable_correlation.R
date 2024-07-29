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
              "Growing degree days (GDD5)", "Minimum winter tempterature (WMIN)",
              "Temperature seasonality (tasCV)", "Soil moisture (soilM)",
              "Total annual rainfall (RAIN)")

# PROCESS DATA --------------------------------------

cov_R_1990 <- lapply(cov_R, function(x) { x[[1]] }) %>% rast(.)
cov_R_2015 <- lapply(cov_R, function(x) { x[[2]] }) %>% rast(.)

names(cov_R_1990) <- names(cov_R_2015) <- covNames

# CORRELATION PLOT -------------------------------------------------

# 1990 correlation plot
cor1990 <- cov_R_1990 %>%
  as.data.frame %>%
  cor
png(filename = '../Writing/Plots/corr_all_1990.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor1990, type = "upper", order = "original", 
         tl.col = "black", addCoef.col="black" ,tl.srt = 45)
dev.off()

# 2015 correlation plot
cor2015 <- cov_R_2015 %>%
  as.data.frame %>%
  cor
png(filename = '../Writing/Plots/corr_all_2015.png',
    width = 30, height = 30, units = "cm", res = 300)
corrplot(cor2015, type = "upper", order = "original" , 
         tl.col = "black",addCoef.col="black" ,tl.srt = 45)
dev.off()

# Combine into single plot
corBoth <- cowplot::ggdraw(clip = "on") +
  cowplot::draw_image('../Writing/Plots/corr_all_1990.png', -0.1, 0.4, 1.25, 0.625) +
  cowplot::draw_image('../Writing/Plots/corr_all_2015.png', -0.1, -0.1, 1.25, 0.625) +
  cowplot::draw_label("(a) 1990", 0.1, 0.95, size = 22) +
  cowplot::draw_label("(b) 2015", 0.1, 0.45, size = 22)

# Save
ggsave(filename = paste0("../Writing/Plots/", "allCorr.png"),
       corBoth,
       dpi = 600,
       units = "px", width = 4000, height = 8400)  

# Remove individual plots as not needed
unlink('../Writing/Plots/corr_all_1990.png')
unlink('../Writing/Plots/corr_all_2015.png') 
