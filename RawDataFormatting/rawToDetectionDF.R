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
# input
community <- c("Sinangoe", "Siona", "Zabalo", "Remolino", "San Pablo")
communityAbrv <- c("SGE", "SNA", "ZAB", "REM", "SPA")
species <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Metachirus nudicaudatus", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Didelphis marsupialis", "Leopardus pardalis")


for (j in 1:length(community)) {
    # clear workspace so that there isn't any accidental overwriting
    rm(list = setdiff(ls(), c("community", "communityAbrv", "species", "j", "i")))

    # load data
    Data <- read.csv(paste0(community[j], "/Data/", communityAbrv[j], "IndependentRecordsFormatted.csv")) # just independent records
    Traps <- read.csv(paste0(community[j], "/Data/", communityAbrv[j], "StationsFormatted.csv"))
    Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))

    # replace all Mazama species with Mazama sp.
    Data$Species <- gsub("Mazama americana", "Mazama sp.", Data$Species)
    Data$Species <- gsub("Mazama nemorivaga", "Mazama sp.", Data$Species)
    Data$Species <- gsub("Mazama gouazoubira", "Mazama sp.", Data$Species)

    # explore total detections
    total <- data.frame(table(Data$Species)) # number of detections / species
    colnames(total) <- c("Species", "Total")

    # get detection history for each species
    for (i in 1:length(species)) {
        # if species is not in the data, skip
        if (length(which(Data$Species == species[i])) == 0) {
            next
        }
        # camera operability matrix
        Operation = cameraOperation(
            CTtable = Traps,
            stationCol = "Station",
            cameraCol = "Camera",
            setupCol = "Setup_date",
            retrievalCol = "Retrieval_date",
            hasProblems = TRUE,
            byCamera = FALSE,
            allCamsOn = FALSE,
            camerasIndependent = FALSE,
            dateFormat = "%Y-%m-%d",
            writecsv = FALSE
        )

        # occasion length
        occasion = 2

        # species detection histories for occupancy analyses
        DetHis = detectionHistory(
            recordTable = Data,
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
            # maxNumberDays = 90, #need to think about this
            species = species[i]
        ) # change species here

        justDetHis <- DetHis[["detection_history"]]

        write.csv(justDetHis,
            paste0(community[j], "/Data/", communityAbrv[j], gsub(" ", "", species[i]), ".csv"),
            row.names = T
        )
    }



}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---------------------- ALL COMMUNITIES TOGETHER -----------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
rm(list = ls())
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
species <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Metachirus nudicaudatus", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Didelphis marsupialis", "Leopardus pardalis")

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

    print(paste0("Finished ", species[i], " of ", length(species), " :)"))
}

