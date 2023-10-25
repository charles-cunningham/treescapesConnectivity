# HEADER --------------------------------------------
#
# Author: Charles Cunningham
# Email: charles.cunningham@york.ac.uk
# 
# Script Name: Run Omniscape 
#
# Script Description: Parallelised run for Omniscape across
# multiple parameters (radius and resistance) and for both 1990
# and 2015.

# SET UP ---------------------------------------------

# N.B. Need to run this line once in Julia before running main script on cluster
# using Pkg; Pkg.add(["Omniscape", "Rasters", "Plots", "DelimitedFiles"])

# Read in packages
using Omniscape, Rasters, Plots, DelimitedFiles

# Each iteration is one row of the input table, specified by job script
i = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

# Read in input table
parameters = readdlm("/users/cac567/scratch/Treescapes/Omniscape/inputTable.txt",
                     Int16, # set table type
                     header = true)

# Get full file path to '.ini' file 'i' based on input table
# N.B. omniscape requires full path
INIpath = string("/users/cac567/scratch/",
                 "/Treescapes/Omniscape/radius",
                 parameters[1][i,1], # input radius
                 "_resistance", 
                 parameters[1][i,3], # input resistance
                 "/omniscapeSettings",
                 parameters[1][i,4], # input year
                 ".ini") 

# Run Omniscape
run_omniscape(string(INIpath))

#####

# N.B. How to set number of threads:

# Linux/Mac -
# export JULIA_NUM_THREADS=4

# Windows -
# set JULIA_NUM_THREADS=4