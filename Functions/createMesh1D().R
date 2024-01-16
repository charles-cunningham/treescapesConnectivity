# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Function Name: createMesh1D
#
# Function Description: Create 1D mesh points for all covariate
# values in a spatRaster for use in inlabru modelling.

# Start function:
# spatRast = any covariate spatRast with 1 to n layers
# nBins = how many bins are there to be in your 1D mesh
createMesh1D <- function(spatRast, nBins) {
  
  # Find min of all values in all spatRast layers
  minRange <- spatRast %>%
    as.vector(.) %>%
    min(., na.rm = TRUE)
  
  # Find max of all values in all spatRast layers
  maxRange<- spatRast %>%
    as.vector(.) %>%
    max(., na.rm = TRUE)
  
  # Create 1D mesh using min and max
  # N.B. this adds on a little to either end of mesh to avoid odd boundary effects
  mesh <- inla.mesh.1d(seq(minRange - (maxRange-minRange)/nBins,
                           maxRange - (maxRange-minRange)/nBins, 
                           length.out = nBins),
                       boundary = "free")
  
  return(mesh)
}