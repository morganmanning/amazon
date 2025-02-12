rm(list = ls())

############################################################################
###################### FORMATTING DETECTION MATRICES #######################
############################################################################
# requirements: 
# notes:


setwd("~/Documents/amazon")

library(openxlsx)
library(camtrapR)
library(dplyr)
require(lubridate)


############################################################################
############################# LOAD DATA ####################################
############################################################################

##### Pick a community
# Sinangoe = SGE
# Siona = SNA
# Siekopai = SKP
# Zabalo = ZAB
# Remolino = REM
# San Pablo = SPA
rm(list = ls())
community <- "Remolino"
communityAbrv <- "REM"
Data <- read.csv(paste0(community, "/Data/", communityAbrv, "IndependentRecordsFormatted.csv")) # just independent records
Traps <- read.csv(paste0(community, "/Data/", communityAbrv, "StationsFormatted.csv")) 
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))

##### Pick species of interest
species <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Metachirus nudicaudatus", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Didelphis marsupialis")
#species <- c("Cuniculus paca", "Mazama americana", "Pecari tajacu", "Psophia crepitans")
# paca = Cuniculus paca
# brocket = Mazama americana
# collared peccary = Dicotyles tajacu 
# trumpeter = Psophia crepitans
# brown four-eyed possum = Metachirus nudicaudatus (#1 species in SGE)
# black agouti = Dasyprocta fuliginosa (#2 species in SGE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ----------------------------- DETECTIONS ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# explore total detections
total <- data.frame(table(Data$Species)) # number of detections / species
colnames(total) <- c("Species", "Total")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- OCCUPANCY SET-UP ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
for (i in 1:length(species)) {
# camera operability matrix
Operation = cameraOperation(CTtable = Traps,
                            stationCol = "Station",
                            cameraCol = "Camera",
                            setupCol = "Setup_date",
                            retrievalCol = "Retrieval_date",
                            hasProblems = TRUE,
                            byCamera = FALSE,
                            allCamsOn = FALSE,
                            camerasIndependent = FALSE,
                            dateFormat = "%Y-%m-%d",
                            writecsv = FALSE)

# occasion length
occasion = 2

# species detection histories for occupancy analyses
DetHis = detectionHistory(recordTable = Data,
                          camOp = Operation,
                          output = "binary", # binary or count
                          stationCol = "Station",
                          speciesCol = "Species",
                          recordDateTimeCol = "DateTimeOriginal",
                          recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                          day1 = "Station",
                          occasionLength = occasion,
                          datesAsOccasionNames = FALSE,
                          timeZone = "America/Guayaquil",
                          includeEffort = TRUE,
                          scaleEffort = FALSE,
                          #maxNumberDays = 90, #need to think about this
                          species = species[i]) #change species here

justDetHis <- DetHis[["detection_history"]]

write.csv(justDetHis, 
          paste0(community, "/Data/", communityAbrv, gsub(" ", "", species[i]), ".csv"), 
          row.names=T)
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---------------------- ALL COMMUNITIES TOGETHER -----------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

Data <- read.csv("Global/Data/AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("Global/Data/AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ----------------------------- DETECTIONS ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# explore total detections
total <- data.frame(table(Data$Species)) # number of detections / species
colnames(total) <- c("Species", "Total")
total <- total[order(-total$Total),]
head(total, 10)

# replace all Mazama species with Mazama sp.
Data$Species <- gsub("Mazama americana", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama nemorivaga", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama gouazoubira", "Mazama sp.", Data$Species)

##### Pick species of interest
species <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Metachirus nudicaudatus", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Didelphis marsupialis")

#species <- c("Cuniculus paca", "Mazama americana", "Pecari tajacu", "Psophia crepitans")
# paca = Cuniculus paca
# brocket = Mazama americana
# collared peccary = Dicotyles tajacu 
# trumpeter = Psophia crepitans
# brown four-eyed possum = Metachirus nudicaudatus (#4 species overall)
# black agouti = Dasyprocta fuliginosa (#1 species overall)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- OCCUPANCY SET-UP ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
for (i in 1:length(species)) {
  # camera operability matrix
  Operation = cameraOperation(CTtable = Traps,
                              stationCol = "Station",
                              cameraCol = "Camera",
                              setupCol = "Setup_date",
                              retrievalCol = "Retrieval_date",
                              hasProblems = TRUE,
                              byCamera = FALSE,
                              allCamsOn = FALSE,
                              camerasIndependent = FALSE,
                              dateFormat = "%Y-%m-%d",
                              writecsv = FALSE)
  
  # occasion length
  occasion = 2
  
  # species detection histories for occupancy analyses
  DetHis = detectionHistory(recordTable = Data,
                            camOp = Operation,
                            output = "binary", # binary or count
                            stationCol = "Station",
                            speciesCol = "Species",
                            recordDateTimeCol = "DateTimeOriginal",
                            recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                            day1 = "Station",
                            occasionLength = occasion,
                            datesAsOccasionNames = FALSE,
                            timeZone = "America/Guayaquil",
                            includeEffort = TRUE,
                            scaleEffort = FALSE,
                            #maxNumberDays = 90, #need to think about this
                            species = species[i]) #change species here
  
  justDetHis <- DetHis[["detection_history"]]
  
  write.csv(justDetHis, 
            paste0("Global/Data/All", gsub(" ", "", species[i]), ".csv"), 
            row.names=T)
}

