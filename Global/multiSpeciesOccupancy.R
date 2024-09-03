##### PAIRWISE CONDITIONAL PROBABILITY COMPARISON #####
setwd("/Users/morganmanning/Documents/amazon/Global/Data")
setwd("~/Documents/amazon/Global/Data")

################################################################################
# ------------------------------ START UP -------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# load necessary packages
require(dplyr)
require(lubridate)
require(camtrapR)
require(unmarked)
require(ggplot2)
require(rphylopic)
require(knitr)
require(kableExtra)

# load in the necessary data
Data <- read.csv("AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
ZABhunting <- read.csv("../../Zabalo/Data/HuntingData2018.csv")

# get tally of each species at each community
speciesTally <- Data |> 
  group_by(Species, CommunityName) |>
  summarize(nDetections = n()) |>
  #filter(nDetections > 10) |>
  group_by(Species) |>
  mutate(nCommunities = n()) # |>
  # filter(nCommunities == 4) # only pull species that were detected in all four communities

huntingTally <- ZABhunting |> 
  group_by(Species) |>
  summarize(nHunted = n()) |>
  filter(nHunted > 10)

################################################################################
# ------------------------- DETECTION MATRICES --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# predator/prey; 
  # ocelot (Leopardus pardalis), uncommon prey (agouti/paca), common prey (green acouchi [Myoprocta pratti], brown four-eyed opossum [Metachirus nudicaudatus], red squirrel [Hadrosciurus spadiceus])
  # include time of day as a covariate (?) for prey since possum is nocturnal, squirrel is diurnal, and ocelot is crepsular

# two most common species (e.g., P(peccary|paca) being lower in disintegrated spaces because people are going to hunt areas with more desirable species; 
  # peccary, paca, black agouti (most spotted but not super hunted, Dasyprocta fuliginosa)

# two species that are totally different with one hunted and one not hunted as a sort of control; 
  # trouble finding two hunted species that aren't predator prey

# two species that are going to be competing in the same niche
  #

# for fun
  # (Panthera onca)

# species of interest                       ************* INPUT ***************
species <- c("Pecari tajacu", 
             "Mazama americana", 
             "Cuniculus paca", 
             "Psophia crepitans",
             "Dasyprocta fuliginosa", # by FAR the most detected species, but not hunted in ZAB much
             "Leopardus pardalis")
commonNames <- c("Collared peccary", 
                 "Red brocket", 
                 "Lowland paca", 
                 "Grey-winged trumpeter", 
                 "Black agouti",
                 "Ocelot") # listTitles

# set up blank lists
detection <- list()

# camera operability matrix
Operation <- cameraOperation(CTtable = Traps,
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

# detection matrices
for (i in 1:length(species)) {
  # occasion length
  occasion = 10 # picked arbitrarily
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
  detection[[i]] <- DetHis[["detection_history"]]
  names(detection)[i] <- species[i]
}










