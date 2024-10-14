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
require(ggplot2)
require(rphylopic)
require(knitr)
require(kableExtra)

# load in the necessary data
functionalTraits <- read.csv("speciesAttributesManualInput.csv")
Data <- read.csv("AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
clusterings <- read.csv("sklearnClusterings.csv")


################################################################################
# ----------------------------- CLUSTERING ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# make the row names = Genus_species
rownames(functionalTraits) <- functionalTraits$Name

# remove M. rufina bc it's not one of the detected species
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
      occasion = 10
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

# prep dataframe for estimate input
estimates <- data.frame()

# run single-species, single-season occupancy model for each cluster
occupancyModelList <- list()
middleman <- list()

for (i in 1:length(allClusterDetHis)){ # take the clusters that are clustered by __
    ## make dataframe to fill in with estimates
    perClustering <- data.frame(clusterBy = rep(clusterRange[i], times = clusterRange[i]),
                        clusterNumber = NA,
                        occupancyEstimate = NA,
                        occupancySE = NA,
                        detectionEstimate = NA,
                        detectionSE = NA)
    
    for (j in 1:length(allClusterDetHis[[i]])){ # take each cluster and run null occupancy on it
        ## run null model for each cluster
        cluster <- allClusterDetHis[[i]][[j]]
        ufo <- unmarkedFrameOccu(cluster,
                                 siteCovs = NULL,
                                 obsCovs = NULL)
        output <- occu(~1 ~1, ufo, 
                       control = 10000,
                       starts = c(0, 0)) # starting values for parameters
        middleman[[j]] <- output

        ## fill in estimates DF
        occPred <- backTransform(output, "state")
        detPred <- backTransform(output, "det")
        #perClustering$clusterBy[j] <- clusterRange[i]
        perClustering$clusterNumber[j] <- j
        perClustering$occupancyEstimate[j] <- occPred@estimate
        perClustering$occupancySE[j] <- SE(occPred)
        perClustering$detectionEstimate[j] <- detPred@estimate
        perClustering$detectionSE[j] <- SE(detPred)
    }
    estimates <- rbind(estimates, perClustering)
    occupancyModelList[[i]] <- middleman # have to do this or get a vector too long error
  print(paste0("Just finished loop ", i, " out of ", length(allClusterDetHis), " :)"))
}



################################################################################
# -------------------------- PLOTTING ESTIMATES -------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# average detection and occupancy per clustering
averagedByClustering <- estimates |>
    group_by(clusterBy) |>
    mutate(averagedOccupancy = mean(occupancyEstimate),
           averagedOccupancySE = mean(occupancySE, na.rm = TRUE),
           averagedDetection = mean(detectionEstimate),
           averagedDetectionSE = mean(detectionSE, na.rm = TRUE)) |>
    mutate(averagedOccupancySE = replace(averagedOccupancySE, is.na(occupancySE), 1),
           averagedDetectionSE = replace(averagedDetectionSE, is.na(detectionSE), 1))

# plotting averaged occupancy estimates
ggplot(averagedByClustering, aes(x = clusterBy,  y = averagedOccupancy)) +
    geom_ribbon(aes(ymin = averagedOccupancy - averagedOccupancySE,
                    ymax = averagedOccupancy + averagedOccupancySE),
                alpha = 0.5, fill = "lightgreen") +
    geom_line(color = "darkgreen") +
    geom_point(color = "forestgreen") +
    coord_cartesian(ylim = (c(0.75,1))) +
    scale_x_discrete(name ="Number of clusters", 
                     limits=seq(min(clusterRange), max(clusterRange), by = 4)) +
    ylab("Averaged occupancy probability estimate") +             
    theme_bw()
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/clustersOccupancy.png", 
       width = 7, height = 5)

# plotting averaged detection estimates
ggplot(averagedByClustering, aes(x = clusterBy, y = averagedDetection)) +
    geom_ribbon(aes(ymin = averagedDetection - averagedDetectionSE,
                    ymax = averagedDetection + averagedDetectionSE),
                alpha = 0.5, fill = "lightblue")  +
    geom_point(color = "steelblue") +
    geom_line(color = "steelblue") +
    coord_cartesian(ylim = (c(0.45,.75))) +
    scale_x_discrete(name ="Number of clusters", 
                     limits=seq(min(clusterRange), max(clusterRange), by = 4)) +
    ylab("Averaged detection probability estimate") +             
    theme_bw()
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/clustersDetection.png", 
       width = 7, height = 5)

# plotting occupancy estimates
ggplot(estimates, aes(x = clusterBy, y = occupancyEstimate)) +
    geom_point() +
    #geom_line() +
    geom_errorbar(aes(ymin = occupancyEstimate - occupancySE,
                    ymax = occupancyEstimate + occupancySE)) +
  ylim(c(0,1)) +
    theme_bw()

# plotting detection estimates
ggplot(estimates, aes(x = as.factor(clusterBy), y = detectionEstimate)) +
    geom_point() +
    #geom_line() +
    geom_errorbar(aes(ymin = detectionEstimate - detectionSE,
                    ymax = detectionEstimate + detectionSE)) +
  ylim(c(0,1)) +
    theme_bw()



################################################################################
# -------------------------- BROCKET CASE STUDY -------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

brocketSpecies <- c("Mazama americana", "Mazama gouazoubira",
                    "Mazama nemorivaga", "Mazama sp.")
brockPic <- get_phylopic(uuid = get_uuid(name = "Mazama americana", n = 1))

# detection matrices for all brockets and all brockets
brocketMatrices <- list()
for (i in 1:length(brocketSpecies)) {
      # occasion length
      occasion = 10
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
                                species = brocketSpecies[i]) #change species here
            brocketMatrices[[i]] <- DetHis[["detection_history"]]
          }
summedMatrices <- Reduce("+", justDetHis)
summedMatrices[summedMatrices > 0] <- 1
brocketMatrices[[5]] <- summedMatrices # matrix for all brockets
names(brocketMatrices) <- c(brocketSpecies, "All brockets")

# run single-species, single-season null occupancy model for each cluster
occuBrocket <- list()
brocketEstimates  <- data.frame(species = c(brocketSpecies, "All brockets"),
                        occupancyEstimate = NA,
                        occupancySE = NA,
                        detectionEstimate = NA,
                        detectionSE = NA)
for (i in 1:length(brocketMatrices)){ # take each cluster and run null occupancy on it
        ## run null model for each cluster
        brocketN <- brocketMatrices[[i]]
        ufo <- unmarkedFrameOccu(brocketN,
                                 siteCovs = NULL,
                                 obsCovs = NULL)
        output <- occu(~1 ~1, ufo, 
                       control = 10000,
                       starts = c(0, 0)) # starting values for parameters
        occuBrocket[[j]] <- output
        ## fill in estimates DF
        occPred <- backTransform(output, "state")
        detPred <- backTransform(output, "det")
        brocketEstimates$occupancyEstimate[i] <- occPred@estimate
        brocketEstimates$occupancySE[i] <- SE(occPred)
        brocketEstimates$detectionEstimate[i] <- detPred@estimate
        brocketEstimates$detectionSE[i] <- SE(detPred)
    }
brocketEstimates

# plot null occupancy estimates
ggplot(brocketEstimates, aes(x = c(as.character(gsub(" ", "\n", brocketEstimates$species))),
                             y = occupancyEstimate)) +
  geom_rect(xmin = 0, xmax = 1.5, ymin = -0.5, ymax = 1.5,
            fill = 'yellow', linetype = "dashed", color = "yellow2", alpha = 0.05) +
  geom_errorbar(aes(ymin = occupancyEstimate - occupancySE,
                    ymax = occupancyEstimate + occupancySE),
                width = 0.2, color  = "darkgreen") +
  add_phylopic(brockPic, alpha = 0.2, x = 4.75, y = 0.85, ysize = 0.3) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen", size = 2) +
  ylim(c(0,1)) +
  ylab("Occupancy probability estimate") +
  xlab("Brocket species") +
  theme_bw()
ggsave("../Figures/brocketClusteredOccupancy.png", width = 7, height = 5)
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/brocketClusteredOccupancy.png", 
       width = 7, height = 5)

# plot null detection estimates
ggplot(brocketEstimates, aes(x = c(as.character(gsub(" ", "\n", brocketEstimates$species))),
                             y = detectionEstimate)) +
  geom_rect(xmin = 0, xmax = 1.5, ymin = -0.5, ymax = 1.5,
            fill = 'yellow',linetype = "dashed", color = "yellow2", alpha = 0.05) +
    geom_errorbar(aes(ymin = detectionEstimate - detectionSE,
                    ymax = detectionEstimate + detectionSE),
               width = 0.2, color  = "dodgerblue") +
    geom_line(color = "dodgerblue") +
    geom_point(color = "dodgerblue", size = 2) +
  add_phylopic(brockPic, alpha = 0.2, x = 4.75, y = 0.85, ysize = 0.3) +
  ylim(c(0,1)) +
    ylab("Detection probability estimate") +
    scale_x_discrete(name = "Brocket species") +
    theme_bw()
ggsave("../Figures/brocketClusteredDetection.png", width = 7, height = 5)
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/brocketClusteredDetection.png", 
       width = 7, height = 5)


# run single-species, single-season occupancy model with covariates
# prediction dataframes
communityPrediction <- data.frame(Community = c("Zábalo", "Siekopai", "Sinangoe", "Siona"),
                                  Year = as.factor(2022))
communityPrediction$Community <- factor(communityPrediction$Community,
                                        levels = c("Zábalo", "Siekopai", "Sinangoe", "Siona"))
yearPrediction <- data.frame(Community = "Siekopai",
                             Year = as.factor(unique(year(Data$DateTimeOriginal))))

# run models
brocketCommPred <- list()
brocketYearPred <- list()
occuCovariatesBrocket <- list()
for (i in 1:length(brocketMatrices)){ # take each cluster and run null occupancy on it
    
        ## run model for each species
    brocketN <- brocketMatrices[[i]]

    ## site covariates
    siteCovariates <- data.frame(Station = rownames(brocketN),
                                 Community = NA,
                                 Year = NA)
    siteCovariates$Community <- ifelse(grepl('^ZAB', siteCovariates$Station), 'Zábalo',
                                ifelse(grepl('^SKP', siteCovariates$Station), 'Siekopai',
                                ifelse(grepl('^SGE', siteCovariates$Station), 'Sinangoe', "Siona")))
    siteCovariates$Community <- factor(siteCovariates$Community,
                                       levels = c("Zábalo", "Siekopai", "Sinangoe", "Siona"))
    siteCovariates$Year <- ifelse(siteCovariates$Community == "Zábalo", "2018", "2022")
    
        ufo <- unmarkedFrameOccu(brocketN,
                                 siteCovs = siteCovariates,
                                 obsCovs = NULL)
        output <- occu(~1 ~Community, ufo, 
                       control = 10000) 
        occuCovariatesBrocket[[j]] <- output
    # occupancy prediction
    brocketCommPred[[i]] <- unmarked::predict(fitList(output), type = "state",
                                              new = communityPrediction, append = TRUE)
    #brocketYearPred[[i]] <- unmarked::predict(fitList(output), type = "state",
                                              #new = yearPrediction, append = TRUE)
    }

names(brocketCommPred) <- names(brocketMatrices)
#names(brocketYearPred) <- names(brocketMatrices)

# make dataframe for community predictions
brocketCommPredDF <- dplyr::bind_rows(brocketCommPred) 
brocketCommPredDF$Species <- rep(names(brocketCommPred), each = nrow(brocketCommPred[[1]]))
brocketCommPredDF

# make dataframe for year predictions
# brocketYearPredDF <- dplyr::bind_rows(brocketYearPred) 
# brocketYearPredDF$Species <- rep(names(brocketYearPred), each = nrow(brocketYearPred[[1]]))
# brocketYearPredDF

# plot occupancy by species and community
colors <- c("Zábalo" = "darkgreen", "Siekopai" = "forestgreen", 
            "Sinangoe" = "yellowgreen", "Siona" = "goldenrod")
ggplot(brocketCommPredDF, aes(x = c(as.character(gsub(" ", "\n", brocketCommPredDF$Species))),
                              y = Predicted, color = Community)) +
  geom_rect(xmin = 0, xmax = 1.5, ymin = -0.5, ymax = 1.5, 
            fill = 'yellow', linetype = "dashed", color = "yellow2", alpha = 0.01) +
  geom_errorbar(aes(ymin = Predicted - SE, ymax = Predicted + SE),
                width = 0.2) +
  geom_point(size = 2) +
  add_phylopic(brockPic, alpha = 0.2, x = 4.75, y = 0.85, ysize = 0.3) +
  ylab("Occupancy probability estimate") +
  scale_x_discrete(name = "Brocket species") +
  scale_color_manual(values = colors) +
  theme_bw()
ggsave("../Figures/brocketCommunityOccupancy.png", width = 7, height = 5)
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/brocketCommunityOccupancy.png", 
       width = 7, height = 5)


## other potentially useful stats for paper
table(Data$CommunityName)



################################################################################
# ----------------------- ABUNDANCE AND DIVERSITY -----------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# remove all unknown species
noUnknownsSGE <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Sinangoe"),]
noUnknownsZAB <- Data[(Data$Species != "NAN NAN" & Data$Species != "NA NA" & Data$CommunityName == "Zabalo"),]
noUnknownsSNA <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Siona"),]
noUnknownsSPA <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "San Pablo"),]
noUnknownsREM <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Remolino"),]


cameraInfo <- Data %>% 
  filter(Species != "N/D N/D" & Species != "NAN NAN" & Species != "NA NA") %>%
  group_by(CommunityName) %>%
  mutate(StartDate = min(DateTimeOriginal),
         EndDate = max(as.Date(DateTimeOriginal)),
         Year = year(DateTimeOriginal)) %>%
  summarise(OperatingDays = round(as.numeric(max(DateTimeOriginal)-
                                               min(DateTimeOriginal))),
            StartDate = min(DateTimeOriginal),
            EndDate = max(DateTimeOriginal),
            numberOfStations = length(unique(Station))) 
            # numberOfCamerasPerStation = round(length(unique(CameraName))/length(unique(Station)))
cameraInfo$StartDate <- format(cameraInfo$StartDate, "%Y-%m-%d")
cameraInfo$EndDate <- format(cameraInfo$EndDate, "%Y-%m-%d")
# put cameraInfo so that it is in the following order, Zabalo, Remolino, Sinangoe, San Pablo, then Siona
cameraInfo <- cameraInfo[c(5,1,3,2,4),] # order by natural area
cameraInfo$CommunityName <- gsub(cameraInfo$CommunityName, pattern = "Zabalo", replacement = "Zábalo")

# per station
  # Sinangoe
siteDiversitySGE <- noUnknownsSGE %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySGE <- siteDiversitySGE %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Zábalo
siteDiversityZAB <- noUnknownsZAB %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversityZAB <- siteDiversityZAB %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Siona
siteDiversitySNA <- noUnknownsSNA %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySNA <- siteDiversitySNA %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Siekopai
siteDiversitySKP <- noUnknownsSKP %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySKP <- siteDiversitySKP %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# per community
wholeDiversitySGE <- noUnknownsSGE %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Sinangoe", PercentNaturalArea = 0.766, 
         OperatingDays = round(as.numeric(max(noUnknownsSGE$DateTimeOriginal)-
                                            min(noUnknownsSGE$DateTimeOriginal)), 3))
wholeDiversityZAB <- noUnknownsZAB %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Zábalo", PercentNaturalArea = 0.936, 
         OperatingDays = round(as.numeric(max(noUnknownsZAB$DateTimeOriginal)-
                                            min(noUnknownsZAB$DateTimeOriginal)), 3))
wholeDiversitySNA <- noUnknownsSNA %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Siona", PercentNaturalArea = 0.754, 
         OperatingDays = round(as.numeric(max(noUnknownsSNA$DateTimeOriginal)-
                                            min(noUnknownsSNA$DateTimeOriginal)), 3))
wholeDiversitySKP <- noUnknownsSKP %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Siekopai", PercentNaturalArea = 0.807, 
         OperatingDays = round(as.numeric(max(noUnknownsSKP$DateTimeOriginal)-
                                      min(noUnknownsSKP$DateTimeOriginal)), 3))

# abundance and diversity for all communities
communityAbundance <- rbind(wholeDiversityZAB, wholeDiversitySKP, wholeDiversitySGE, wholeDiversitySNA)

communityDiversity <- communityAbundance %>%
  group_by(Community, PercentNaturalArea) %>%
  summarise(nIndiv=sum(abundance),
            nSpecies = length(unique(Species)),
            OperatingDays = mean(OperatingDays),
            shannonIndex = round(-sum((abundance/sum(abundance))*log(abundance/sum(abundance))), 3),
            simpsonIndex = round(1-sum((abundance/sum(abundance))^2), 3)) 
communityDiversity <- arrange(communityDiversity, desc(PercentNaturalArea))
communityDiversity$Community <- factor(communityDiversity$Community, 
                                       levels = communityDiversity$Community)
communityDiversity


ggplot(communityDiversity, aes(x = Community, y = PercentNaturalArea, fill = Community)) +
  geom_bar(stat="identity") +
  ylab("Percent natural surrounding area") +
  scale_fill_manual(values = colors) +
  ylim(c(0,1)) +
  theme_bw()
ggsave("../Figures/percentNatArea.png", width = 7, height = 5)
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/percentNatArea.png", 
       width = 7, height = 5)


################################################################################
# ----------------------------- MAKE TABLES -----------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

head(communityDiversity)
kbl(communityDiversity, col.names = c("Community", "Percent Natural Area", 
                                      "Number of Detections", "Number of Species",
                                      "Number of Sampling Days",
                                      "Shannon Diversity Index", "Simpson Diversity Index")) %>%
  kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
  save_kable(file = "../Figures/communityDiversity.png", zoom = 1.5)

kbl(communityDiversity[,c("Community", "PercentNaturalArea", "shannonIndex", "simpsonIndex")], 
    col.names = c("Community", "Percent Natural Area", 
                  "Shannon Diversity Index", "Simpson Diversity Index")) %>%
  kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
  save_kable(file = "../Figures/communityDiversitySummary.png", zoom = 2)

kbl(communityDiversity, col.names = c("Community", "Percent Natural Area", 
                                      "Number of Detections", "Number of Species",
                                      "Number of Sampling Days",
                                      "Shannon Diversity Index", "Simpson Diversity Index")) %>%
  kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
  save_kable(file = "~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/communityDiversity.png",
             zoom = 1.5)

## Camera info
head(cameraInfo)
kbl(cameraInfo, col.names = c("Community", "Number of Sampling Days", 
                              "Sampling Start Date", "Sampling End Date",
                              "Number of Sites")) %>%
  kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
  save_kable(file = "../Figures/siteInfo.png",
             zoom = 1.5)
kbl(cameraInfo, col.names = c("Community", "Number of Sampling Days", 
                                      "Sampling Start Date", "Sampling End Date",
                                      "Number of Sites",
                                      "Number of Cameras per Site")) %>%
  kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
  save_kable(file = "~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/siteInfo.png",
             zoom = 1.5)
 



################################################################################
# ----------------------- BROCKET DETECTION PLOTS -----------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

brocketN <- brocketMatrices[['All brockets']]
ufo <- unmarkedFrameOccu(brocketN,
                         siteCovs = NULL,
                         obsCovs = NULL)
colnames(ufo@y) <- 1:ncol(ufo@y)
rownames(ufo@y) <- 1:nrow(ufo@y)
require(reshape2)
meltedComb <- melt(ufo@y)
meltedComb$value <- as.factor(meltedComb$value)

ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
  geom_tile(colour = "gray50") +
  scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
  scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  xlab(paste("Time (1 unit = ~", 
             10*2,
             # number of columns/time steps divided by the clumping factor, times 2
             "days)")) +
  ylab("Camera trap site") +
  #ggtitle("Grouped brocket detection matrix") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size = 25, hjust = 0.5))
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/groupedBrocketDetection.png", 
       width = 5, height = 5)

# one species of brocket
brocketN <- brocketMatrices[['Mazama americana']]
ufo <- unmarkedFrameOccu(brocketN,
                         siteCovs = NULL,
                         obsCovs = NULL)
colnames(ufo@y) <- 1:ncol(ufo@y)
rownames(ufo@y) <- 1:nrow(ufo@y)
require(reshape2)
meltedComb <- melt(ufo@y)
meltedComb$value <- as.factor(meltedComb$value)

ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
  geom_tile(colour = "gray50") +
  scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
  scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  xlab(paste("Time (1 unit = ~", 
             10*2,
             # number of columns/time steps divided by the clumping factor, times 2
             "days)")) +
  ylab("Camera trap site") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(plot.title = element_text(size = 25, hjust = 0.5))
ggsave("~/Dropbox/UF/Spring2024/WIS6505CQuantitativeAnalysis/Final/MamericanaDetection.png", 
       width = 5, height = 5)

