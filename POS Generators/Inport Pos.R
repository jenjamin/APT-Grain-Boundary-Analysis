#### Inputs .pos file of GB region that has been determined by the analysis ####

setwd(dirname(parent.frame(2)$ofile))

readposR <- function(posFileName) {
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  fileInfo <- file.info(posFileName)
  fileSize <- fileInfo$size / sizeOfFloat
  to.read = file(posFileName, "rb")
  posFile <- readBin(to.read,
                     double(),
                     size = 4,
                     n = fileSize,
                     endian = "big")
  close(to.read)
  posFile <-
    t(matrix(posFile, nrow = 4, dimnames = list(c("x", "y", "z", "m"))))
  posFile <- as.data.frame(posFile)
  return(posFile)
}

PosFile <- readposR("GB Containing ROI.pos")

#### Input RangeFile and Tidy ####
RangeFile <- read.delim("Test Range.rrng", colClasses = "character", header = FALSE)
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

#### Assigning ion based on m ####
IonList <- lapply(PosFile$m, function(MassToCharge) RangesDF2$Ion[between(MassToCharge,
                                                                          RangesDF2$Start,
                                                                          RangesDF2$End)])
IonList <- lapply(IonList, function(x) if(identical(x, character(0))) NA_character_ else x)
PosFile$Ion <- unlist(IonList)
rm(IonList)
PosFile$Ion[is.na(PosFile$Ion)] = "Unranged"

Ion <- data.frame(matrix(NA,
                         nrow = nrow(PosFile)))
for(i in seq(1,nrow(RangesDF2),1)){
  Name <- RangesDF2$Ion[i]
  Ion <<- cbind(
    Ion,
    PosFile %>%
      mutate(Name = ifelse(PosFile$m > RangesDF2$Start[i] &
                             PosFile$m < RangesDF2$End[i],
                           RangesDF2$Ion[i], NA)) %>%
      select(Name))
}
Ion$Noise <- "Noise"
PosFile$Ion <- apply(Ion, 1, function(x) na.omit(x)[1])
rm(i, Ion)

#### Plot Mass Spec ####

myTheme <- function(relSize = 16) {
  return(theme_classic(relSize))
}

# Create colour list from range file
myColors <- c(paste0("#",
                     (RangesDF2 %>%
                        select(Color, Ion) %>%
                        distinct())$Color),
              "#000000")
names(myColors) <- c((RangesDF2 %>%
                        select(Color, Ion) %>%
                        distinct())$Ion,
                     "Unranged")

ggplot(PosFile %>%
         filter(m < 33 & m > 0)) +
  geom_histogram(aes(m, fill = Ion),
                 binwidth = 0.01) + 
  scale_y_continuous(trans = log10_trans()) +
  scale_fill_manual(values = myColors) +
  myTheme() +
  labs(x = "Mass-To-Charge-State Ratio (m/z)",
       y = "Count") +
  theme(legend.position = "bottom")
