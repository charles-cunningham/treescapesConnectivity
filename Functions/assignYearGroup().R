# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Function Name: assignYearGroup
#
# Function Description: Group a numeric vector (years)
# into discrete time periods. Must be 2 groups.

# Create function to assign group to species records
assignYearGroup <- function(value, 
                            yearGroupRange1 = c(NA,NA),
                            yearGroupRange2 = c(NA,NA)) {
  
  if (value >= yearGroupRange1[1] & value < yearGroupRange1[2] & !is.na(value)) {
    return(1)
  }
  
  else if (value >= yearGroupRange2[1] & value < yearGroupRange2[2] & !is.na(value)) {
    return(2)
  }
  else
    (return(NA))
}