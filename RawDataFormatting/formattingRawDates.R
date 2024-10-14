setwd("~/Documents/amazon/RawDataFormatting")
require(lubridate)
require(dplyr)
##### Cleaning up dates, etc. on raw datasets 

### Zabalo
rm(list=ls())
Data <- read.csv("../Zabalo/Data/ZABIndependentRecords.csv") # just independent records
Traps <- read.csv("../Zabalo/Data/ZABStations.csv") 

Data$Species <- gsub("-", " ", Data$Species)
Data$Species <- gsub("  ", " ", Data$Species)
Data$Species <- gsub("Dicotyles tajacu", "Pecari tajacu", Data$Species)
Data$Species <- gsub("Hadroscirus spadiceus", "Hadrosciurus spadiceus", Data$Species)
Data$Species <- gsub("Tinamu major", "Tinamus major", Data$Species)
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, "%m/%d/%y %H:%M")
Data$Date <- as.Date(Data$Date, format = "%m/%d/%y")
Data$DateTimeOriginal <- with(Data, ymd(Date) + hms(Time))
Data$Station <- paste0("ZAB", Data$Station)
Data$Camera <- paste0("ZAB", Data$Camera)

Traps$Setup_date <-as.Date(Traps$Setup_date, format = "%m/%d/%y")
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, format = "%m/%d/%y")
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, "%m/%d/%y")
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, "%m/%d/%y")
Traps$Problem2_from <- parse_date_time(Traps$Problem2_from, "%m/%d/%y")
Traps$Problem2_to <- parse_date_time(Traps$Problem2_to, "%m/%d/%y")
Traps$Problem3_from <- parse_date_time(Traps$Problem3_from, "%m/%d/%y")
Traps$Problem3_to <- parse_date_time(Traps$Problem3_to, "%m/%d/%y")
Traps <- subset(Traps, Setup_date != Retrieval_date) # need this so don't get error in conversion to dectection matrix
Traps$Station <- paste0("ZAB", Traps$Station)
Traps$Camera <- paste0("ZAB", Traps$Camera)

Data$CommunityName <- "Zabalo"
names(Data)[names(Data) == 'Camera'] <- 'CameraName'
names(Traps)[names(Traps) == 'x'] <- 'gps_x'
names(Traps)[names(Traps) == 'y'] <- 'gps_y'
Traps$CommunityName <- "Zabalo"

write.csv(Data, "../Zabalo/Data/ZABIndependentRecordsFormatted.csv")
write.csv(Traps, "../Zabalo/Data/ZABStationsFormatted.csv")

### Sinangoe
rm(list=ls())
Data <- read.csv("../Sinangoe/Data/SGEIndependentRecords.csv") # just independent records
Traps <- read.csv("../Sinangoe/Data/SGEStations.csv") 

Data$Species <- gsub("  ", " ", Data$Species)
Data$Species <- gsub("Dicotyles tajacu", "Pecari tajacu", Data$Species)
Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))

Data$CommunityName <- "Sinangoe"
Traps$CommunityName <- "Sinangoe"

write.csv(Data, "../Sinangoe/Data/SGEIndependentRecordsFormatted.csv")
write.csv(Traps, "../Sinangoe/Data/SGEStationsFormatted.csv")


### Siona
rm(list=ls())
Data <- read.csv("../Siona/Data/SNAIndependentRecords.csv") # just independent records
Traps <- read.csv("../Siona/Data/SNAStations.csv") 

Data$Species <- gsub("  ", " ", Data$Species)
Data$Species <- gsub("Dicotyles tajacu", "Pecari tajacu", Data$Species)
Traps$Setup_date <-as.Date(Traps$Setup_date, "%d/%m/%Y")
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, "%d/%m/%Y")
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, "%d/%m/%Y")
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, "%d/%m/%Y")
Traps$Problem2_from <- parse_date_time(Traps$Problem2_from, "%d/%m/%Y")
Traps$Problem2_to <- parse_date_time(Traps$Problem2_to, "%d/%m/%Y")
Traps$Camera <- paste0("SNA", Traps$Camera)

Data$CommunityName <- "Siona"
Traps$CommunityName <- "Siona"
names(Traps)[names(Traps) == 'X'] <- 'Obs'

write.csv(Data, "../Siona/Data/SNAIndependentRecordsFormatted.csv")
write.csv(Traps, "../Siona/Data/SNAStationsFormatted.csv")


### Siekopai
rm(list=ls())
Data <- read.csv("../Siekopai/Data/SKPIndependentRecords.csv") # just independent records
Traps <- read.csv("../Siekopai/Data/SKPStations.csv")

Data$Species <- gsub("  ", " ", Data$Species)
Data$Species <- gsub("Dicotyles tajacu", "Pecari tajacu", Data$Species)
Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Camera <- paste0("SKP", Traps$Camera)

Data$CommunityName <- "Siekopai"
Traps$CommunityName <- "Siekopai"
names(Traps)[names(Traps) == 'X'] <- 'Obs'

write.csv(Data, "../Siekopai/Data/SKPIndependentRecordsFormatted.csv")
write.csv(Traps, "../Siekopai/Data/SKPStationsFormatted.csv")


### All records
rm(list=ls())
SNArecords <- read.csv("../Siona/Data/SNAIndependentRecordsFormatted.csv") # just independent records
ZABrecords <- read.csv("../Zabalo/Data/ZABIndependentRecordsFormatted.csv") # just independent records
SGErecords <- read.csv("../Sinangoe/Data/SGEIndependentRecordsFormatted.csv") # just independent records
SKPrecords <- read.csv("../Siekopai/Data/SKPIndependentRecordsFormatted.csv") # just independent records

SNAstations <- read.csv("../Siona/Data/SNAStationsFormatted.csv")
ZABstations <- read.csv("../Zabalo/Data/ZABStationsFormatted.csv")
SGEstations <- read.csv("../Sinangoe/Data/SGEStationsFormatted.csv")
SKPstations <- read.csv("../Siekopai/Data/SKPStationsFormatted.csv")

ZABrecords$Station <- as.character(ZABrecords$Station)
ZABstations$Station <- as.character(ZABstations$Station)
SNAstations$Camera <- as.character(SNAstations$Camera)
SKPstations$Camera <- as.character(SKPstations$Camera)


# make master dataframes
allCommunityRecords <- dplyr::bind_rows(SNArecords,
                                        ZABrecords,
                                        SGErecords,
                                        SKPrecords)
allCommunityRecords <- allCommunityRecords %>%
  select(Station, CameraName, Species, DateTimeOriginal, CommunityName)

allCommunityStations <- dplyr::bind_rows(SNAstations,
                                        ZABstations,
                                        SGEstations,
                                        SKPstations)
allCommunityStations <- allCommunityStations %>%
  select(CommunityName, Station, Camera, gps_y, gps_x, Setup_date, Retrieval_date, 
         Problem1_from, Problem1_to, Problem2_from, Problem2_to, Problem3_from, Problem3_to,
         Total, Obs)

# need to divide Siekopai into Remolino and San Pablo
# after looking at communities on QGIS and the territory shapefiles/vectors provided by Bob in the Teams drive, need to reclassify some communities
# San Pablo is SNA3 and then SKP31-37
allCommunityRecords$CommunityName[allCommunityRecords$Station == "SNA3"] <- "San Pablo"
allCommunityStations$CommunityName[allCommunityStations$Station == "SNA3"] <- "San Pablo"
for(i in 31:37){
  site_to_change <- paste0("SKP", i)
  allCommunityRecords$CommunityName[allCommunityRecords$Station == site_to_change] <- "San Pablo"
  allCommunityStations$CommunityName[allCommunityStations$Station == site_to_change] <- "San Pablo"
}

# Remolino is SKP1-30
for(i in 1:30){
  site_to_change <- paste0("SKP", i)
  allCommunityRecords$CommunityName[allCommunityRecords$Station == site_to_change] <- "Remolino"
  allCommunityStations$CommunityName[allCommunityStations$Station == site_to_change] <- "Remolino"
}

write.csv(allCommunityRecords, "../Global/Data/AllIndependentRecordsFormatted.csv")
write.csv(allCommunityStations, "../Global/Data/AllStationsFormatted.csv")


###### Subset San Pablo and Remolino to have independent detection histories
# load in the community data
allCommunityRecords <- read.csv("../Global/Data/AllIndependentRecordsFormatted.csv")
allCommunityRecords$X <- NULL
allCommunityStations <- read.csv("../Global/Data/AllStationsFormatted.csv")
allCommunityStations$X <- NULL
allCommunityStations$Obs <- NULL

# subset
SPArecords <- allCommunityRecords[allCommunityRecords$CommunityName == "San Pablo",]
REMrecords <- allCommunityRecords[allCommunityRecords$CommunityName == "Remolino",]
SPAstations <- allCommunityStations[allCommunityStations$CommunityName == "San Pablo",]
REMstations <- allCommunityStations[allCommunityStations$CommunityName == "Remolino",]

# remove unnecessary rows (blank columns affection creation of cameraOperation())
REMstations$Problem2_from <- NULL
REMstations$Problem2_to <- NULL
REMstations$Problem3_from <- NULL
REMstations$Problem3_to <- NULL
SPAstations$Problem2_from <- NULL
SPAstations$Problem2_to <- NULL
SPAstations$Problem3_from <- NULL
SPAstations$Problem3_to <- NULL

# save it
write.csv(SPArecords, "../San Pablo/Data/SPAIndependentRecordsFormatted.csv")
write.csv(SPAstations, "../San Pablo/Data/SPAStationsFormatted.csv")

write.csv(REMrecords, "../Remolino/Data/REMIndependentRecordsFormatted.csv")
write.csv(REMstations, "../Remolino/Data/REMStationsFormatted.csv")




