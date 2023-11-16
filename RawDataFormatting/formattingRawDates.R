setwd("~/Documents/amazon/RawDataFormatting")

##### Cleaning up dates, etc. on raw datasets 

### Zabalo
Data <- read.csv("Zabalo/ZABIndependentRecords.csv") # just independent records
Traps <- read.csv("Zabalo/ZABStations.csv") 

Data$Species <- gsub("-", " ", Data$Species)
Data$Species <- gsub("Pecari tajacu", "Dicotyles tajacu", Data$Species)
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
Traps <- subset(Traps, Setup_date != Retrieval_date)

write.csv(Data, "../Zabalo/Data/ZABIndependentRecordsFormatted.csv")
write.csv(Data, "Zabalo/ZABIndependentRecordsFormatted.csv")

write.csv(Traps, "Zabalo/ZABStationsFormatted.csv")
write.csv(Traps, "../Zabalo/Data/ZABStationsFormatted.csv")

### Sinangoe
Data <- read.csv("Sinangoe/SGEIndependentRecords.csv") # just independent records
Traps <- read.csv("Sinangoe/SGEStations.csv") 

Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))

write.csv(Data, "Sinangoe/SGEIndependentRecordsFormatted.csv")
write.csv(Data, "../Sinangoe/Data/SGEIndependentRecordsFormatted.csv")

write.csv(Traps, "Sinangoe/SGEStationsFormatted.csv")
write.csv(Traps, "../Sinangoe/Data/SGEStationsFormatted.csv")


### Siona
Data <- read.csv("Siona/SNAIndependentRecords.csv") # just independent records
Traps <- read.csv("Siona/SNAStations.csv") 

Traps$Setup_date <-as.Date(Traps$Setup_date, "%d/%m/%Y")
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, "%d/%m/%Y")
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, "%d/%m/%Y")
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, "%d/%m/%Y")
Traps$Problem2_from <- parse_date_time(Traps$Problem2_from, "%d/%m/%Y")
Traps$Problem2_to <- parse_date_time(Traps$Problem2_to, "%d/%m/%Y")

write.csv(Data, "Siona/SNAIndependentRecordsFormatted.csv")
write.csv(Data, "../Siona/Data/SNAIndependentRecordsFormatted.csv")

write.csv(Traps, "Siona/SNAStationsFormatted.csv")
write.csv(Traps, "../Siona/Data/SNAStationsFormatted.csv")



### Siekopai
Data <- read.csv("Siekopai/SKPIndependentRecords.csv") # just independent records
Traps <- read.csv("Siekopai/SKPStations.csv")

Traps$Setup_date <-as.Date(Traps$Setup_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Retrieval_date <- as.Date(Traps$Retrieval_date, tryFormats = c("%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"))
Traps$Problem1_from <- parse_date_time(Traps$Problem1_from, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))
Traps$Problem1_to <- parse_date_time(Traps$Problem1_to, c("%d/%m/%Y", "%d/%m/%y", "%m/%d/%Y", "%m/%d/%y"))

write.csv(Data, "Siekopai/SKPIndependentRecordsFormatted.csv")
write.csv(Data, "../Siekopai/Data/SKPIndependentRecordsFormatted.csv")

write.csv(Traps, "Siekopai/SKPStationsFormatted.csv")
write.csv(Traps, "../Siekopai/Data/SKPStationsFormatted.csv")






