# Script for analysing a grain boundary containing pos file
# Will determine start and end of interface, GB composition and width
# Will return a .pos file of just the GB region

require(tidyverse)

#### User defined  variables - change these for each analysis ####
# Path to .pos file
PosLocation <- c("Test Pos and Range Files/Simulated GB.pos")
# Path to .rrng file
RangeFileLocation <- c("Test Pos and Range Files/Range For Simulation.rrng")
# Detection Efficiency of instrument
DetectionEfficiency = 0.37
# Orientation is the "x", "y", or "z" direction that is perpendicular to the interface
Orientation = "z"
# Number of bins you would like the dataset to be divided into
Bins = 40

#### Set WD to source file location ####
setwd(dirname(parent.frame(2)$ofile))

#### Create new directory for generated files ####
NewDirectory <- file.path(dirname(getwd()), basename(getwd()), paste0("Output Files - ",gsub(":","",Sys.time())))
dir.create(NewDirectory)

#### Read intial .pos file ####
source("Scripts For Grain Boundary Analysis/readposR.R")
PosFileInput <- readposR(PosLocation)

##### Combine .pos file with .rrng file to get chemical information ####
source("Scripts For Grain Boundary Analysis/Combine Pos and Range File.R")
PosFileRanger(PosFileInput,
              RangeFileLocation)

rm(Ion, PosFileInput, PosFileRanger, RangeFileLocation, readposR)

#### Generate 1D concentration profile ####
source("Scripts For Grain Boundary Analysis/PosOneDimensionalPlot.R")
OneDCountFunc(PosFile = RangedPos,
              NumberOfBins = Bins,
              Direction = Orientation)

# Calculate extent of GB region
source("Scripts For Grain Boundary Analysis/GBExtentCalculator.R")

# Export .pos of GB region
source("Scripts For Grain Boundary Analysis/GB Pos Generator.R")

# Contour plots etc. of GB region 
