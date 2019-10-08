#Calculates Gibbs Excess of Interface

GibbsExcessCalculator <- function(Area, DetectionEfficiency) {
  DividingSurfacePosition <- mean(GB$Distance) / max(GrainB$Distance)
  LowerSurfacePosition <- min(GB$Distance) / max(GrainB$Distance)
  UpperSurfacePosition <- max(GB$Distance) / max(GrainB$Distance)
  
  TotalVolumeConc <- as.data.frame(colSums(
    ElementCountDF %>%
      select(-Distance, -Atom_Count)
    )) %>%
    rownames_to_column("Element") %>%
    rename(NumberCounts = 2)
  
  TotalNIons <- sum(TotalVolumeConc$NumberCounts)
  
  TotalVolumeConc <- TotalVolumeConc %>%
    mutate(GlobalFraction = NumberCounts / TotalNIons)
  
  GrainAConc <- as.data.frame(colSums(GrainA %>%
                                        select(-Distance, -Atom_Count))) %>%
    rownames_to_column("Element") %>%
    rename(NumberCounts = 2)
  
  GrainAIons <- sum(GrainAConc$NumberCounts)
  
  GrainAConc <- GrainAConc %>%
    transmute(GrainAFraction = NumberCounts / GrainAIons)
  
  GrainBConc <- as.data.frame(colSums(GrainB %>%
                                        select(-Distance, -Atom_Count))) %>%
    rownames_to_column("Element") %>%
    rename(NumberCounts = 2)
  
  GrainBIons <- sum(GrainBConc$NumberCounts)
  
  GrainBConc <- GrainBConc %>%
    transmute(GrainBFraction = NumberCounts / GrainBIons)
  
  
  GibbsExcess <<- (cbind(TotalVolumeConc,
                         GrainAConc,
                         GrainBConc)) %>%
    mutate(
      GrainAGE = (1 / (Area * DetectionEfficiency)) * TotalNIons * (GlobalFraction - GrainAFraction),
      GrainBGE = (1 / (Area * DetectionEfficiency)) * TotalNIons * (GlobalFraction - GrainBFraction),
      GibbsExcess = (1 / (Area * DetectionEfficiency)) * TotalNIons * (
        GlobalFraction -
          GrainAFraction * (DividingSurfacePosition) -
          GrainBFraction * (1 - DividingSurfacePosition)
      ),
      GibbsExcessLower = (1 / (Area * DetectionEfficiency)) * TotalNIons * (
        GlobalFraction -
          GrainAFraction * (LowerSurfacePosition) -
          GrainBFraction * (1 - LowerSurfacePosition)
      ),
      GibbsExcessUpper = (1 / (Area * DetectionEfficiency)) * TotalNIons * (
        GlobalFraction -
          GrainAFraction * (UpperSurfacePosition) -
          GrainBFraction * (1 - UpperSurfacePosition)
      )
      
    )
  
}
