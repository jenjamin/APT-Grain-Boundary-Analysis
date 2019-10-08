# Script for creating simulated pos file of data
library(tidyverse)

# Set WD to source file location
setwd(dirname(parent.frame(2)$ofile))

POSGenerator <- function(AtomicDensity, XRange, YRange, ZRange,
                         CompositionTable){
  
  TotalNumberAtoms = AtomicDensity * 
    (max(XRange)-min(XRange)) * 
    (max(YRange)-min(YRange)) *
    (max(ZRange)-min(ZRange))
  x <- sample(seq(min(XRange), max(XRange), 0.00001), TotalNumberAtoms, replace = FALSE)
  y <- sample(seq(min(YRange), max(YRange), 0.00001), TotalNumberAtoms, replace = FALSE)
  z <- sample(seq(min(ZRange), max(ZRange), 0.00001), TotalNumberAtoms, replace = FALSE)
  
  SimulatedPos <- data.frame(x,y,z) %>%
    mutate(m = sample(CompositionTable$Mass, n(), prob = CompositionTable$Abundance, replace = TRUE))
  
}

#### Generate matrix pos ####

MatrixA <- POSGenerator(AtomicDensity = 40, 
                        XRange = c(-10, 10),
                        YRange = c(-10, 10), 
                        ZRange = c(-15, -1),
                        CompositionTable = data.frame("Element" = c("Fe", "Ni", "P2Si3"),
                                                      "Abundance"= c(98.9,1.0,0.1),
                                                      "Mass" = c(1, 2, 3)))

MatrixB <- POSGenerator(AtomicDensity = 40, 
                        XRange = c(-10, 10),
                        YRange = c(-10, 10), 
                        ZRange = c(1, 15),
                        CompositionTable = data.frame("Element" = c("Fe", "Ni", "P2Si3"),
                                                      "Abundance"= c(98.9,1.0,0.1),
                                                      "Mass" = c(1, 2, 3)))
#### GB pos generation ####
GB <- POSGenerator(AtomicDensity = 40, 
                   XRange = c(-10, 10),
                   YRange = c(-10, 10), 
                   ZRange = c(-1, 1),
                   CompositionTable = data.frame("Element" = c("Fe", "Ni", "P2Si3"),
                                                 "Abundance"= c(90.0,6.0,4.0),
                                                 "Mass" = c(1, 2, 3)))

#### Create overall pos ####
SimulatedPos <- rbind(MatrixA, MatrixB, GB)
#Add gaussian noise to z position
SimulatedPos <- SimulatedPos %>%
  mutate(z = z + rnorm(n(),0,0.5)) %>%
  filter(-10 < z & z < 10)

rm(MatrixA, MatrixB, GB)

source("../Scripts For Grain Boundary Analysis/writeposR.R")

writeposR(SimulatedPos, "../Test Pos and Range Files//Simulated GB.pos")
