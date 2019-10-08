# Create function to generate pos file
# Input atomic denisty
# ROI dimensions
library(tidyverse)

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

#### Input RangeFile and Tidy ####
RangeFile <- read.delim("../Test Pos and Range Files/Range For Simulation.rrng", colClasses = "character", header = FALSE)
RowsToSkip <- as.numeric(gsub("Number=", "", RangeFile[2,])) + 5
Ranges <- RangeFile %>% slice(RowsToSkip:n())
rm(RowsToSkip, RangeFile)

RangesDF <- data.frame()

i = 0
for(i in unique(str_count(Ranges$V1, ":") - 2)){
  
  Elements <- c()
  for(j in seq(1,i,1)){
    Elements <- append(Elements,paste("Element",j))
  }
  
  ColumnNames <- c("Start", "End", "Volume",
                   Elements, "Color")
  
  Ranges %>%
    mutate(NumberIons = str_count(V1, ":") - 2) %>%
    filter(NumberIons == i) %>%
    separate(V1,
             ColumnNames,
             sep = " ")
  
  RangesDF <- bind_rows(RangesDF, 
                        Ranges %>%
                          mutate(NumberIons = str_count(V1, ":") - 2) %>%
                          filter(NumberIons == i) %>%
                          separate(V1,
                                   ColumnNames,
                                   sep = " ")
  )
  
}

rm(ColumnNames, Elements,i, j) 

#### Creat R-friendly range file ####

RangesDF2 <- cbind(
  RangesDF %>%
    mutate(Start = as.numeric(str_extract(Start,"[^=]+$")),
           End = as.numeric(str_extract(End, "[^=]+$")),
           Volume = gsub("Vol:", "", Volume),
           Color = gsub("Color:", "", Color)) %>%
    select(Start, End, Volume, Color),
  RangesDF %>%
    select(contains("Element")) %>%
    unite("Ion") %>%
    mutate(Ion = paste(gsub("1|Name|:|NA|_| ","",Ion)))
)

rm(Ranges, RangesDF)

#### Range pos####
Ion <- data.frame(matrix(NA,
                          nrow = nrow(SimulatedPos)))
for(i in seq(1,nrow(RangesDF2),1)){
  Name <- RangesDF2$Ion[i]
  Ion <<- cbind(
    Ion,
    SimulatedPos %>%
    mutate(Name = ifelse(SimulatedPos$m > RangesDF2$Start[i] &
                        SimulatedPos$m < RangesDF2$End[i],
                        RangesDF2$Ion[i], NA)) %>%
    select(Name))
}
Ion$Noise <- "Noise"
SimulatedPos$Ion <- apply(Ion, 1, function(x) na.omit(x)[1])
rm(i, Ion)


#### Plot one D conc plot ####
ROILength = max(SimulatedPos$z) - min(SimulatedPos$z)

OneDConcPlot <- SimulatedPos %>% 
  group_by(Distance = cut(z, breaks= seq(floor(min(SimulatedPos$z)),
                                         ceiling(max(SimulatedPos$z)),
                                         by = 0.1)),
           Ion) %>%
  summarise(Ioncount = n()) %>%
  ungroup() %>%
  spread(Ion, Ioncount) %>%
  mutate(Distance = as.numeric(as.character((ceiling(ROILength))*(row_number()/n()))))

OneDConcPlot[is.na(OneDConcPlot)] <- 0

ggplot(OneDConcPlot) +
  geom_point(aes(Distance, P2Si3))

source("../Scripts For Grain Boundary Analysis/writeposR.R")

writeposR(SimulatedPos %>% select(-Ion), "Simulated GB.pos")
