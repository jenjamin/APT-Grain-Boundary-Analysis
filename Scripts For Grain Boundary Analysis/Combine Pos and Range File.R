# Function to combine .pos file and .rrng file

PosFileRanger <- function(PosFile, RangeFilePath){
  #### Input RangeFile and Tidy ####
  RangeFile <- read.delim(RangeFilePath, colClasses = "character", header = FALSE)
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
  
  RangeInfo <<- cbind(
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
  RangeInfo$Ion <- gsub(",","", RangeInfo$Ion)
  
  rm(Ranges, RangesDF)
  
  #### Range pos####
  Ion <<- as.vector(matrix(NA,
                           nrow = nrow(PosFile)))
  
  for (IonName in 1:length(RangeInfo$Ion)) {
    IsIonName <- with(PosFile,
                      m >= slice(RangeInfo, IonName)$Start &
                        m <= slice(RangeInfo, IonName)$End)
    IonLocations <- which(IsIonName == TRUE)
    Ion <- replace(Ion, IonLocations, slice(RangeInfo, IonName)$Ion)
  }
  
  # Replace NA with "Noise"
  Ion <- replace(Ion, which(is.na(Ion) == TRUE), "Noise")
  
  PosFile["Ion"] <- Ion
  assign("RangedPos", PosFile, envir=.GlobalEnv)
  
  rm(Ion)
}




