# Script for plotting variation of concentration or count with distance in .pos file
myTheme <- function(relSize = 16) {
  return(theme_classic(relSize))
}

OneDCountFunc <- function(PosFile, NumberOfBins, Direction){
  NoiselessPos <- PosFile %>%
    filter(Ion != "Noise")
  
  QuantiledDF <- NoiselessPos %>%
    mutate(
      Quant =  ntile(NoiselessPos[,Direction], NumberOfBins),
      DirectionDistance = NoiselessPos[,Direction]
      )
  
  OneDCount <<- merge(
    QuantiledDF %>%
      group_by(Quant) %>%
      summarise(Distance = mean(DirectionDistance)) %>%
      ungroup(),
    QuantiledDF %>%
      group_by(Quant, Ion) %>%
      summarise(Ioncount = n()) %>%
      ungroup() %>%
      spread(Ion, Ioncount)) %>%
    select(-Quant)
  
  OneDCount[is.na(OneDCount)] <<- 0
  
  OneDConc <<- OneDCount %>%
    gather(Ion, Count, -Distance) %>%
    group_by(Distance) %>%
    mutate(Total = sum(Count),
           Concentration = 100 * Count/Total) %>%
    select(-Count, -Total) %>%
    spread(Ion, Concentration)
  
  # Create colour list from range file
  myColors <- c(paste0("#",
                       (RangeInfo %>%
                          select(Color, Ion) %>%
                          distinct())$Color),
                "#000000")
  names(myColors) <- c((RangeInfo %>%
                          select(Color, Ion) %>%
                          distinct())$Ion,
                       "Unranged")
  
  print(ggplot(OneDCount %>%
           gather(Ion, Count, - Distance)) +
    geom_point(aes(Distance, Count, color = Ion)) +
    scale_color_manual(values = myColors) +
    myTheme() +
    labs(x = "Distance (nm)",
         y = "Count") +
    theme(legend.position = "bottom"))
  
  print(ggplot(OneDConc %>%
                 gather(Ion, Concentration, - Distance)) +
          geom_point(aes(Distance, Concentration, color = Ion)) +
          scale_color_manual(values = myColors) +
          myTheme() +
          labs(x = "Distance (nm)",
               y = "Concentration (Ionic %)") +
          theme(legend.position = "bottom"))
  

    
  }

