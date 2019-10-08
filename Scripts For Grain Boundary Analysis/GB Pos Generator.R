# Script to select grain boundary region from pos and then export .pos file of region

GrainBoundaryPoS <- RangedPos %>%
  select(-Ion) %>%
  filter(UQ(as.symbol(Orientation)) >= InterfaceStart &
           UQ(as.symbol(Orientation)) <= InterfaceEnd)

source("Scripts For Grain Boundary Analysis/writeposR.R")

setwd(NewDirectory)
writeposR(GrainBoundaryPoS,
          "Grain Boundary Region.pos")

write.table(
  paste0("These files were generated from a pos file at the following location - ",
         suppressWarnings(normalizePath(PosLocation)), "\n",
         "The width of the interface is: ",round(GBWidth,2)," nm (2dp) \n",
         "The area used for Gibbs Excess calculations was ", Area," nm^2"),
  file = "Grain Boundary Analysis Summary.txt", sep = "\t",
            row.names = FALSE)
setwd('..')
