
# Goal: make plots to show the number of detections for each species and how often they overlap at stations

setwd("/Users/morganmanning/Documents/amazon/Global/Data")
setwd("~/Documents/amazon/Global/Data")

################################################################################
################################################################################
################################################################################

# load necessary packages
require(dplyr)
library(tidyverse)
require(lubridate)
require(camtrapR)
require(unmarked)
require(ggplot2)
require(knitr)
require(kableExtra)
require(stringr)
require(gridExtra)
require(ggpubr)
require(reshape2)

# load in the necessary data
Data <- read.csv("AllIndependentRecordsFormatted.csv")
Traps <- read.csv("AllStationsFormatted.csv")
covariates <- read.csv("AllCommunityCovariates.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))

################################################################################
################################################################################
################################################################################

# replace all Mazama species with Mazama sp.
Data$Species <- gsub("Mazama americana", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama nemorivaga", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama gouazoubira", "Mazama sp.", Data$Species)
Data <- Data[!Data$Species %in% c("N/D N/D", "NA NA", "NAN NAN"),]

# unique species
species <- unique(Data$Species)

### GOAL: make a dataframe with the number stations each species were detected at, and the number of stations that both species were detected at
# make a dataframe to hold the pairwise data
detections <- as.data.frame(t(combn(species, 2)))
colnames(detections) <- c("Species1", "Species2")

# make columns for the number of sites each species was detected at 
detections$Species1_Sites <- NA
detections$Species2_Sites <- NA
detections$Both_Species_Sites <- NA
detections$Species1_Sites <- sapply(detections$Species1, function(x) length(unique(Data$Station[Data$Species == x])))
detections$Species2_Sites <- sapply(detections$Species2, function(x) length(unique(Data$Station[Data$Species == x])))

# add a column for the number of sites that both species overlapped at
for (i in 1:nrow(detections)) {
  # get the number of sites that both species were detected at
  sp1Sites <- unique(Data$Station[Data$Species == detections$Species1[i]])
  sp2Sites <- unique(Data$Station[Data$Species == detections$Species2[i]])
  nOverlap <- length(intersect(sp1Sites, sp2Sites))
  detections$Both_Species_Sites[i] <- nOverlap
}

# make a column for each station for each species that has the number of times each species was detected
detections[, paste(unique(Data$Station), "Species1", sep = "_")] <- NA
detections[, paste(unique(Data$Station), "Species2", sep = "_")] <- NA

# loop through each station and each species and fill in the dataframe
for (i in 1:nrow(detections)) { # takes ~5 minutes
    for (j in 1:length(unique(Data$Station))) {
        # get the number of detections for each species at each station
        detections[i, paste(unique(Data$Station)[j], "Species1", sep = "_")] <- nrow(Data[Data$Species == detections$Species1[i] & Data$Station == unique(Data$Station)[j], ])
        detections[i, paste(unique(Data$Station)[j], "Species2", sep = "_")] <- nrow(Data[Data$Species == detections$Species2[i] & Data$Station == unique(Data$Station)[j], ])
    }

    # print a message every 100 iterations
    if (i %% 100 == 0) {
        print(paste("Finished", i, "of", nrow(detections), "rows :)"))
    }
}


# find the number of times each species was detected within the same time step using detection matrices for each species
# camera operability matrix
Operation <- cameraOperation(
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

# make a detection matrix for each species
detHistList <- list()
for (i in 1:length(species)) {
    # occasion length
    occasion <- 2 # picked arbitrarily
    # species detection histories for occupancy analyses
    DetHis <- detectionHistory(
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
    detHistList[[i]] <- DetHis[["detection_history"]]
    names(detHistList)[i] <- species[i]
}


# make a dataframe that has the number of stations each species was detected at, and then a column for the number of times species A and species B were detected at the same station within the same time step
timeSpecificDetections <- as.data.frame(t(combn(species, 2)))
colnames(timeSpecificDetections) <- c("Species1", "Species2")

# make columns to fill in
timeSpecificDetections$Species1_nDet <- NA
timeSpecificDetections$Species2_nDet <- NA
timeSpecificDetections$Species1_nSites <- NA
timeSpecificDetections$Species2_nSites <- NA
timeSpecificDetections$Species1_nTimeSteps <- NA
timeSpecificDetections$Species2_nTimeSteps <- NA
timeSpecificDetections$nOverlapTimeSteps <- NA

for (i in 1:nrow(timeSpecificDetections)){
    # get the total number of raw detections for each species
    timeSpecificDetections$Species1_nDet[i] <- sum(Data$Species == timeSpecificDetections$Species1[i])
    timeSpecificDetections$Species2_nDet[i] <- sum(Data$Species == timeSpecificDetections$Species2[i])
    
    # get the number of stations each species was detected at
    timeSpecificDetections$Species1_nSites[i] <- length(unique(Data$Station[Data$Species == timeSpecificDetections$Species1[i]]))
    timeSpecificDetections$Species2_nSites[i] <- length(unique(Data$Station[Data$Species == timeSpecificDetections$Species2[i]]))

    # get the total number of time steps for each species
    timeSpecificDetections$Species1_nTimeSteps[i] <- sum(detHistList[[timeSpecificDetections$Species1[i]]], na.rm = TRUE)
    timeSpecificDetections$Species2_nTimeSteps[i] <- sum(detHistList[[timeSpecificDetections$Species2[i]]], na.rm = TRUE)

    # get the number of time steps that both species overlapped
    timeSpecificDetections$nOverlapTimeSteps[i] <- sum(detHistList[[timeSpecificDetections$Species1[i]]] & detHistList[[timeSpecificDetections$Species2[i]]], na.rm = TRUE)

    # print a message every 100 iterations
    if (i %% 100 == 0) {
        print(paste("Finished", i, "of", nrow(timeSpecificDetections), "rows :)"))
    }

}

# make a matrix that has the number of times each species was detected at the same time step
pairwiseDetectionMatrix <- matrix(nrow = length(species), ncol = length(species), data = NA, 
                                   dimnames = list(species, species)) 

for (i in 1:nrow(timeSpecificDetections)){
    # get the number of times each species was detected at the same time step
    sp1 <- timeSpecificDetections$Species1[i]
    sp2 <- timeSpecificDetections$Species2[i]
    nOverlap <- timeSpecificDetections$nOverlapTimeSteps[i]
    pairwiseDetectionMatrix[sp1, sp2] <- nOverlap
    pairwiseDetectionMatrix[sp2, sp1] <- nOverlap

}

## convert to tibble, add row identifier, and shape "long"
pairwiseDetectionMatrix2 <-
    pairwiseDetectionMatrix %>%
    as_tibble() %>%
    mutate(Species1 = rownames(pairwiseDetectionMatrix)) %>%
    pivot_longer(-Species1, names_to = "Species2", values_to = "value") %>%
    arrange(desc(value)) %>%
    mutate(
        Species1 = as.factor(Species1),
        Species2 = as.factor(Species2)
    ) 

# plot the pairwise detection matrix
ggplot(pairwiseDetectionMatrix2, aes(Species1, Species2)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = value)) +
    scale_fill_gradient(low = "white", high = "red") +
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
    file = paste0("../Figures/MultispeciesModeling/pairwiseDetectionMatrix.png"),
    width = 14, height = 14
)

# remove reciprocal rows
d1 <- as.data.frame(t(apply(timeSpecificDetections, 1, sort)))
timeSpecificDetections_noDups <- timeSpecificDetections[!duplicated(d1), ] %>%
    arrange(desc(nOverlapTimeSteps))



# make a table of the top pairwise detections
kbl(head(timeSpecificDetections_noDups[,c("Species1", "Species2", "Species1_nTimeSteps","Species2_nTimeSteps", "nOverlapTimeSteps")], 30)) %>%
    kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
    save_kable(file = "../Figures/MultispeciesModeling/topPairwiseDetectionsTable.png", zoom = 2)

# species with the most overlapped detections
speciesOfInterest <- c("Mazama sp.", "Dasyprocta fuliginosa", "Cuniculus paca", "Pecari tajacu", "Dasypus novemcinctus", "Psophia crepitans", "Tinamus major", "Metacirus nudicaudatus", "Didelphis marsupialis", "Leopardus pardalis")

# only pull out the rows that contain the species of interest
speciesOfInterestDF <- timeSpecificDetections_noDups[timeSpecificDetections_noDups$Species1 %in% speciesOfInterest & timeSpecificDetections_noDups$Species2 %in% speciesOfInterest, ]
head(speciesOfInterestDF)


# make a faceted plot using timeSpecificDetections_noDups with bars for each species' number of timestep detections and also the number of overlapping detections
# put the number of detections on the y axis and then three bars side by side per faceted plot: 1 bar for species 1 time step detections, 1 bar for species 2 time step detections, and 1 bar for the number of overlapping detections
forPlotting <- melt(speciesOfInterestDF[, c("Species1", "Species2", "Species1_nTimeSteps", "Species2_nTimeSteps", "nOverlapTimeSteps")])
forPlotting <- forPlotting %>%
    mutate(ComboSpecies = paste(Species1, Species2, sep = "_"))
ggplot(forPlotting, aes(x = variable)) +
    facet_wrap(~ComboSpecies) +
    geom_bar()

ggplot(forPlotting, aes(x = variable, y = value)) +
    facet_wrap(~ComboSpecies) +
    geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Species", y = "Number of detections (time step = 2 days)") +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) +
    scale_fill_manual(labels = c("Species 1", "Species 2", "Overlap"), values = c("Species1_nTimeSteps" = "#e3e34b", "Species2_nTimeSteps" = "red", "nOverlapTimeSteps" = "orange")) +
    geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.25) 

ggsave(
    file = "../Figures/MultispeciesModeling/pairwiseDetectionBarPlot.png",
    width = 14, height = 14
)

