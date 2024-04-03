setwd("~/Documents/amazon/RawDataFormatting")
require(lubridate)
require(dplyr)
##### Cleaning up dates, etc. on raw datasets 

### Zabalo
rm(list=ls())
Data <- read.csv("../Zabalo/Data/ZABIndependentRecords.csv") # just independent records
Traps <- read.csv("../Zabalo/Data/ZABStations.csv") 

Data$Species <- gsub("-", " ", Data$Species)
Data$Species <- gsub("Pecari tajacu", "Dicotyles tajacu", Data$Species)
Data$Species <- gsub("Hadroscirus spadiceus", "Hadrosciurus spadiceus", Data$Species)
Data$Species <- gsub("Tinamu major", "Tinamus major", Data$Species)
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, "%m/%d/%y %H:%M")
Data$Date <- as.Date(Data$Date, format = "%m/%d/%y")
Data$DateTimeOriginal <- with(Data, ymd(Date) + hms(Time))
Traps$Setup_date <-as.Date(Traps$Setup_date, format = "%m/%d/%y")
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, format = "%m/%d/%y")
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, "%m/%d/%y")
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, "%m/%d/%y")
Traps$Problem2_from <- parse_date_time(Traps$Problem2_from, "%m/%d/%y")
Traps$Problem2_to <- parse_date_time(Traps$Problem2_to, "%m/%d/%y")
Traps$Problem3_from <- parse_date_time(Traps$Problem3_from, "%m/%d/%y")
Traps$Problem3_to <- parse_date_time(Traps$Problem3_to, "%m/%d/%y")
# Traps <- subset(Traps, Setup_date != Retrieval_date)

Data$Community <- "Zabalo"
names(Data)[names(Data) == 'Camera'] <- 'CameraName'
names(Traps)[names(Traps) == 'x'] <- 'gps_x'
names(Traps)[names(Traps) == 'y'] <- 'gps_y'
Traps$Community <- "Zabalo"


write.csv(Data, "../Zabalo/Data/ZABIndependentRecordsFormatted.csv")
write.csv(Traps, "../Zabalo/Data/ZABStationsFormatted.csv")

### Sinangoe
rm(list=ls())
Data <- read.csv("../Sinangoe/Data/SGEIndependentRecords.csv") # just independent records
Traps <- read.csv("../Sinangoe/Data/SGEStations.csv") 

Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))

Data$Community <- "Sinangoe"
Traps$Community <- "Sinangoe"

write.csv(Data, "../Sinangoe/Data/SGEIndependentRecordsFormatted.csv")
write.csv(Traps, "../Sinangoe/Data/SGEStationsFormatted.csv")


### Siona
rm(list=ls())
Data <- read.csv("../Siona/Data/SNAIndependentRecords.csv") # just independent records
Traps <- read.csv("../Siona/Data/SNAStations.csv") 

Traps$Setup_date <-as.Date(Traps$Setup_date, "%d/%m/%Y")
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, "%d/%m/%Y")
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, "%d/%m/%Y")
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, "%d/%m/%Y")
Traps$Problem2_from <- parse_date_time(Traps$Problem2_from, "%d/%m/%Y")
Traps$Problem2_to <- parse_date_time(Traps$Problem2_to, "%d/%m/%Y")

Data$Community <- "Siona"
Traps$Community <- "Siona"
names(Traps)[names(Traps) == 'X'] <- 'Obs'

write.csv(Data, "../Siona/Data/SNAIndependentRecordsFormatted.csv")
write.csv(Traps, "../Siona/Data/SNAStationsFormatted.csv")


### Siekopai
rm(list=ls())
Data <- read.csv("../Siekopai/Data/SKPIndependentRecords.csv") # just independent records
Traps <- read.csv("../Siekopai/Data/SKPStations.csv")

Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))

Data$Community <- "Siekopai"
Traps$Community <- "Siekopai"
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
  select(Station, CameraName, Species, DateTimeOriginal, Community)

allCommunityStations <- dplyr::bind_rows(SNAstations,
                                        ZABstations,
                                        SGEstations,
                                        SKPstations)
allCommunityStations <- allCommunityStations %>%
  select(Community, Station, Camera, gps_y, gps_x, Setup_date, Retrieval_date, 
         Problem1_from, Problem1_to, Problem2_from, Problem2_to, Problem3_from, Problem3_to,
         Total, Obs)

write.csv(allCommunityRecords, "../Global/Data/allCommunityRecordsFormatted.csv")
write.csv(allCommunityStations, "../Global/Data/allCommunityStationsFormatted.csv")


