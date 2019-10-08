# Function for determining start at end of the interface 

require(tidyverse)
require(PeriodicTable)
require(InterpretMSSpectrum)

#### Creating dataframe after decomposing ions ####
IonList <- colnames(OneDCount %>%
                      select(-Distance))

CorrectDF <- data.frame()
for(Ions in IonList){
  SplitIon <- CountChemicalElements(Ions)
  IonCount <- OneDCount %>%
    select(Distance, Ions)
  CorrectElementCountInIon = data.frame()
  for (Element in 1:length(SplitIon)) {
    CorrectElementCount <- cbind(IonCount %>% select(Distance),
                                 ((IonCount %>%
                                     select(-Distance))* SplitIon[Element]))
    colnames(CorrectElementCount) <- c("Distance", names(SplitIon[Element]))
    CorrectElementCount <- CorrectElementCount %>%
      gather(Element, Count, -Distance)
    CorrectElementCountInIon <<- rbind(CorrectElementCountInIon,CorrectElementCount)
  }
  CorrectDF <<- rbind(CorrectDF, 
                      CorrectElementCountInIon)
}

ElementCountDF <- CorrectDF %>%
  group_by(Element, Distance) %>%
  summarise(Count = sum(Count)) %>%
  spread(Element, Count) %>%
  ungroup() 

ElementCountDF <- ElementCountDF %>%
  mutate(Atom_Count = rowSums(.[2:ncol(ElementCountDF)]))

rm(IonList, Element, Ions, SplitIon, CorrectDF,
   CorrectElementCount, CorrectElementCountInIon) # Remove vars from global envir

print(ggplot(ElementCountDF %>%
               gather(Element, value, -Distance,-Atom_Count) %>%
               filter(Element != "O") %>%
               mutate(`Concentration (at.%)` = 100 * (value/Atom_Count)),
             aes(Distance, `Concentration (at.%)`)) +
        myTheme() + 
        geom_line(aes(colour = Element))
)

#### Calculating Start and End of GB ####

DeterminingElementForThresholdCalculations <- function(){
  IndividualPlot<-readline(prompt="Would you like to plot an individual element (Y/N)? " )
  IndividualPlot<<-as.character(IndividualPlot)
  while(IndividualPlot == c("Y")){
    ElementName<-readline(prompt="Enter element name you would like to plot: " )
    ElementName<<-as.character(ElementName)
    
    if(ElementName %in% colnames(ElementCountDF) == TRUE){
      print(paste0("Plotting Graph for ",ElementName))
    }else{
      print(paste0(ElementName, " is not present in your data. Please select an element that is from the list: ",
                   paste(as.character(colnames(ElementCountDF %>% select(-"Distance", -"Atom_Count"))),collapse= ", "))) 
    }
    
    print(ggplot(ElementCountDF %>% 
                   gather(Element, value, -Distance,-Atom_Count) %>% 
                   filter(Element == ElementName),
                 aes(Distance, value)) +
            geom_line(colour = "black") +
            ylab("Count") +
            ggtitle(paste0(ElementName, " vs Distance")) +
            myTheme())
    IndividualPlot<-readline(prompt="Would you like to plot an individual element (Y/N)? " )
    IndividualPlot<<-as.character(IndividualPlot)
  }
  
  ElementName<-readline(prompt="Enter element you would like to use to calculate thresholds: " )
  ElementName<<-as.character(ElementName)
  PValue<<-as.numeric(readline(prompt="Enter p-value you would like to use to calculate thresholds (typically 0.05 or 0.01): " ))
  
  while(PValue <= 0 | PValue > 1){
    print(paste0("This is an invalid p-value."))
    PValue<<-as.numeric(readline(prompt="Enter p-value you would like to use to calculate thresholds (typically 0.05 or 0.01): " ))
  }
  
  if(ElementName %in% colnames(ElementCountDF) == TRUE){
    print(paste0(ElementName, " will be used to calculate thresholds with a p-value of ",PValue))
  }else{
    print(paste0(ElementName, " is not present in your data")) 
    stop()
  }
}

ThresholdSuitableCheck <-function(){
  ThresholdApprove<-readline(prompt="Does this element give sensible threshold? (Y/N): " )
  ThresholdApprove<<-as.character(ThresholdApprove)
  if(ThresholdApprove == c("Y")){
    
  }else{
    DeterminingElementForThresholdCalculations()
    ThresholdCheckFunction()
  }
}

ThresholdCheckFunction <- function(){
  GrainADistance<<-as.numeric(readline(prompt="Distance up to which the Grain/Phase A is definitely present: " ))
  GrainBDistance<<-as.numeric(readline(prompt="Distance after which the Grain/Phase B is definitely present: " ))
  
  AdjustedPValue <- 1 - (1 - PValue)^(1/nrow(ElementCountDF))
  x <- 1
  Prob <- 1
  while(Prob >= AdjustedPValue){
    Prob <- ppois(x, lambda = mean((ElementCountDF %>% filter(Distance <= GrainADistance))[[ElementName]]), lower=FALSE)
    x <- x + 1
  }
  CutOffAValue <<- x
  
  y <- 1
  Prob <- 1
  while(Prob >= AdjustedPValue){
    Prob <- ppois(y, lambda = mean((ElementCountDF %>% filter(Distance >= GrainBDistance))[[ElementName]]), lower=FALSE)
    y <- y + 1
  }
  CutOffBValue <<- y
  
  print(ggplot(ElementCountDF %>%
                 gather(Element, `Counts in Bin`, -Distance,-Atom_Count) %>%
                 filter(Element == ElementName),
               aes(Distance, `Counts in Bin`)) +
          geom_line(colour = "black") +
          geom_hline(yintercept = CutOffAValue,
                     colour="blue",
                     linetype="dashed") +
          geom_hline(yintercept = CutOffBValue,
                     colour="red",
                     linetype="dashed") +
          geom_vline(xintercept = min((ElementCountDF %>%
                                         filter(UQ(as.symbol(ElementName)) >= CutOffAValue))$Distance),
                     colour="blue") +
          geom_vline(xintercept = max((ElementCountDF %>%
                            filter(UQ(as.symbol(ElementName)) >= CutOffBValue))$Distance),
                     colour="red") +
          geom_text(aes(min((ElementCountDF %>%
                               filter(UQ(as.symbol(ElementName)) >= CutOffAValue))$Distance)/2,
                        CutOffAValue,label="Phase A cutoff value", vjust = -0.5),
                    colour="blue") +
          geom_text(aes((max((ElementCountDF)$"Distance")-
                           (max((ElementCountDF)$"Distance")-
                              max((ElementCountDF %>%
                                     filter(UQ(as.symbol(ElementName)) >= CutOffBValue))$Distance))/2),
                        CutOffBValue,label="Phase B cutoff value", vjust = -0.5),
                    colour="red") +
          geom_text(aes(min((ElementCountDF %>%
                               filter(UQ(as.symbol(ElementName)) >= CutOffAValue))$Distance),
                        0,label="Interface Start",hjust=1.2),
                    colour="blue") +
          geom_text(aes(max((ElementCountDF %>%
                               filter(UQ(as.symbol(ElementName)) >= CutOffBValue))$Distance),
                        0,label="Interface End",hjust=-0.2),
                    colour="red") +
          myTheme() +
          ggtitle(paste0(ElementName," Counts vs Distance")) +
          theme(plot.title = element_text(hjust = 0.5))
  )
  
  ThresholdSuitableCheck()
  
}

DeterminingElementForThresholdCalculations()

ThresholdCheckFunction()

#### Defining Start and End of GB ####
GrainA <- ElementCountDF %>% filter(Distance < min((ElementCountDF %>%
                                                      filter(UQ(as.symbol(ElementName)) >= CutOffAValue))$Distance))
GrainB <- ElementCountDF %>% filter(Distance > max((ElementCountDF %>%
                                                      filter(UQ(as.symbol(ElementName)) >= CutOffBValue))$Distance))
GB <- ElementCountDF %>% 
  filter(Distance < min(GrainB$Distance)) %>% 
  filter(Distance > max(GrainA$Distance))

Area <- (max((RangedPos %>% select(-Orientation, -m, -Ion))[1]) -
           min((RangedPos %>% select(-Orientation, -m, -Ion))[1]) ) * 
  (max((RangedPos %>% select(-Orientation, -m, -Ion))[2]) -
     min((RangedPos %>% select(-Orientation, -m, -Ion))[2]))
source("Scripts For Grain Boundary Analysis/GE Calculator.R")
GibbsExcessCalculator(Area,DetectionEfficiency) #ROI Area and machine detection efficiency


GrainBoundarySumm <- function() {
  GBWidth <<- max(GB$Distance) - min(GB$Distance)
  
  GBComp <- as.data.frame(colSums(GB %>%
                                    select(-Distance,-Atom_Count)))
  
  Total = sum(GBComp)
  
  GBComp <- GBComp %>%
    rownames_to_column("Element") %>%
    mutate(Composition = 100 * `colSums(GB %>% select(-Distance, -Atom_Count))` /
             Total)
  
  GBSummary <<- cbind(GBComp,
                      GibbsExcess %>%
                        select(-Element)) %>%
    select(Element,
           Composition,
           GibbsExcess,
           GibbsExcessLower,
           GibbsExcessUpper,
           GrainAFraction,
           GrainBFraction) %>%
    mutate(
      LowerErr = abs(GibbsExcess - GibbsExcessLower),
      UpperErr = abs(GibbsExcess - GibbsExcessUpper)
    ) %>%
    transmute(
      Element = Element,
      `Interface Composition (at.%)` = round(Composition,4),
      `Gibbs Excess` = GibbsExcess,
      `Gibbs Error` = pmax(LowerErr, UpperErr),
      `Phase A Composition (at.%)` = round(100 * GrainAFraction,4),
      `Phase B Composition (at.%)` = round(100 * GrainBFraction,4)
    )
  
  print(GBSummary)
  setwd(NewDirectory)
  write.csv(GBSummary,"Grain Boundary Calculations.csv",
            row.names = FALSE)
  setwd('..')
  print(paste0("Width of interface is: ",round(GBWidth,2)," nm (2dp)"))
  print(paste0("Area used to calculated interface values was : ",round(Area,2)," nm^2 (2dp)"))
}
GrainBoundarySumm()

InterfaceStart <- min((ElementCountDF %>%
                         filter(UQ(as.symbol(ElementName)) >= CutOffAValue))$Distance)

InterfaceEnd <- max((ElementCountDF %>%
                       filter(UQ(as.symbol(ElementName)) >= CutOffBValue))$Distance)

rm(CutOffAValue, CutOffBValue, ElementCountDF, ElementName, GB, GibbsExcess,
   GibbsExcessCalculator, GrainA, GrainADistance, GrainB, GrainBDistance, IndividualPlot,
   IonCount, PValue, ThresholdApprove, ThresholdCheckFunction, ThresholdSuitableCheck)
# rm(list = c(objects())[!(c(objects()) %in% c("InterfaceEnd","InterfaceStart","GBSummary"))])
