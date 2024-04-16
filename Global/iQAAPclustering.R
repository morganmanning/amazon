##### CALCULATE FUNCTIONAL DIVERSITY #####
setwd("/Users/morganmanning/Documents/amazon/Global/Data")
setwd("~/Documents/amazon/Global/Data")

################################################################################
# ------------------------------ START UP -------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# load necessary packages
require(FD)
require(dplyr)
require(cluster)
require(lubridate)
require(camtrapR)
require(unmarked)

# load in the necessary data
functionalTraits <- read.csv("speciesAttributesManualInput.csv")
Data <- read.csv("AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))



################################################################################
# ----------------------------- CLUSTERING ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# make the row names = Genus_species
rownames(functionalTraits) <- functionalTraits$Name

# remove M. rufina
functionalTraits <- subset(functionalTraits, Name != "Mazama_rufina")

# only select the traits data you're interested in comparing between species
traits <- functionalTraits %>%
  select(Genus,
         Order, 
         Family, 
         ActivityPeriod, 
         BodySizeG, 
         Locomotion, 
         FeedingHabit, 
         FeedingTechnique) 
str(traits)

# calculate functional distance (characters will convert to factors with gowdis())
FDist <- as.matrix(gowdis(traits))

# pull the most similar 5 animals for each animal and put it in a df
nSimilar <- 5
mostSimilar <- data.frame()
for (i in 1:ncol(FDist)){
  speciesOfInterest <- data.frame(species = colnames(FDist)[i],
                                  nearest = rownames(FDist),
                                  functionalDist = FDist[,i])
  speciesOfInterest <- speciesOfInterest[order(speciesOfInterest$functionalDist, 
                                               decreasing = FALSE),] 
  speciesOfInterestTop <- speciesOfInterest[2:(nSimilar+1),] # don't pull comparison of same spp.
  mostSimilar <- rbind(mostSimilar, speciesOfInterestTop)
}
rownames(mostSimilar) <- 1:nrow(mostSimilar)

# k-means clustering
nClusters <- 5
pam <- pam(x = FDist, k = nClusters, diss = TRUE)
clusters <- pam$clustering
clusterDF <- as.data.frame(clusters)
clusterDF <- data.frame(species = gsub("_", " ", rownames(clusterDF)),
                        cluster = clusterDF$clusters)


################################################################################
# ------------------------- DETECTION MATRICES --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# set up blank lists
justDetHis <- list()
clusterDetHis <- list()
allClusterDetHis <- list()

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

# detection matrices for each cluster for each number of clusters
clusterRange <- 3:(nrow(functionalTraits)-1) # 3 because three classes (mammals, birds, reptile)
# Number of clusters 'k' must be in {1,2, .., n-1}; hence n >= 2

for (k in 1:length(clusterRange)) {
  # k-means clustering
  nClusters <- clusterRange[k]
  pam <- pam(x = FDist, k = nClusters, diss = TRUE)
  clusters <- pam$clustering
  clusterDF <- as.data.frame(clusters)
  clusterDF <- data.frame(species = gsub("_", " ", rownames(clusterDF)),
                          cluster = clusterDF$clusters)
  
  for (j in 1:nClusters) { # make detection matrices
    # pull each cluster
    clusterSpecies <- subset(clusterDF, cluster == j)[,1]
    
    for (i in 1:length(clusterSpecies)) {
      
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
                                species = clusterSpecies[i]) #change species here
      
      justDetHis[[i]] <- DetHis[["detection_history"]]
      
    }
    
    summedMatrices <- Reduce("+", justDetHis)
    summedMatrices[summedMatrices > 0] <- 1
    clusterDetHis[[j]] <- summedMatrices # matrices for each cluster
    
  }
  
  allClusterDetHis[[k]] <- clusterDetHis # matrices for each cluster at each grouping factor
  print(paste0("Just finished clustering by ", clusterRange[k], " out of ", max(clusterRange), " :)"))
  
}



################################################################################
# --------------------------- OCCUPANCY MODELS --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# run single-species, single-season occupancy model for each cluster

occupancyModelList <- as.list(rep(0, length(allClusterDetHis)))

for (i in 1:length(allClusterDetHis)){ # take the clusters that are clustered by __
  for (j in 1:length(allClusterDetHis[[i]])){ # take each cluster and run null occupancy on it
    ufo <- unmarkedFrameOccu(allClusterDetHis[[i]][[j]],
                             siteCovs = NULL,
                             obsCovs = NULL)
    occupancyModelList[[i]][[j]] <- occu(~1 ~1, ufo, 
                                         control = 10000,
                                         starts = c(0, 0)) # starting values for parameters
  }
  print(paste0("Just finished loop ", i, " out of ", length(allClusterDetHis), " :)"))
}

## TESTING
ufo <- unmarkedFrameOccu(allClusterDetHis[[54]][[1]],
                         siteCovs = NULL,
                         obsCovs = NULL)
occupancyModelList[[54]][[1]] <- occu(~1 ~1, ufo, 
                                     control = 10000,
                                     starts = c(0, 0))
dim(allClusterDetHis[[54]][[1]]) # 120x66
# ERROR: long vectors not supported yet: /Volumes/Builds/R4/R-4.3.2/src/main/subassign.c:1841







# TO DO:
  # X include all of above in for loop to run all of it for each number of clusters
  # - run occupancy on each cluster for each cluster delineation
  # - plot occupancy predictions for each of said occupancy models





 