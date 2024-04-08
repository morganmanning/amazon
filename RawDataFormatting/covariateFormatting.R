
######################## SETTING UP SITE COVARIATES ##########################
# set my working directory
# setwd("~/Documents/amazon/Zabalo/Data")

# load required packages
require(ggplot2)
require(dplyr)
require(sf)
require(mapview)
require(geosphere)

# load files
siteCovariate <- read.csv("siteCovs2018Raw.csv")
#stations <- read.csv("Stations2018Raw.csv")
stations <- read.csv("ZABStationsFormatted.csv")
stations <- stations %>% 
  distinct(gps_x, gps_y, .keep_all = TRUE)
globalRecords <- read.csv("../../Global/Data/AllStationsFormatted.csv")
trapRecords <- read.csv("ZABIndependentRecordsFormatted.csv") # trapped species records
hunting <- read.csv("HuntingData2018.csv") # harvest species, location, and date
ZABBuffer <- read.csv("../../Global/Data/ZABBuffer.csv", header = FALSE)
SGEBuffer <- read.csv("../../Global/Data/SGEBuffer.csv", header = FALSE)
SNABuffer <- read.csv("../../Global/Data/SNABuffer.csv", header = FALSE)
SKPBuffer <- read.csv("../../Global/Data/SKPBuffer.csv", header = FALSE)



################################################################################
######################## FORMAT HUNTING INTENSITY ##############################
################################################################################

# formatting the data as spatial objects
# latitude = y, longitude = x
cameras <- st_as_sf(stations, coords = c("gps_x", "gps_y"),  crs = 4326) # EPSG:24878 = UTM zone 18S
harvests <- st_as_sf(hunting, coords = c("x", "y"),  crs = 4326)
mapview(harvests, col.regions = "red", alpha = 0.1) +
  mapview(cameras, alpha = 1)

# just the coordinates
cameraCoords <- stations[,c('gps_x', 'gps_y')]
huntingCoords <- hunting[,c('x', 'y')]

# distance from each trap to each hunting occasion
distMat <- distm(as.matrix(huntingCoords), as.matrix(cameraCoords), fun = distGeo)

# weight each hunting occasion for each camera trap (closer hunts weighted more)
expDistances <- matrix(nrow = nrow(hunting), ncol = length(unique(stations$Station)))
for (i in 1:length(unique(stations$Station))) {
  for (j in 1:nrow(hunting)) {
    expDistances[j,i] <- exp(-(distMat[j,i]/1000)) # weight each hunt
    # hunts closer to the camera (i) will be weighted higher
    # all are essentially 0 when not divided by 1000
    # you idiot, you're just turning the distances from meters to kilometers, so it's fine
  }
}
huntingIntensity <- data.frame(station = stations$Station,
                               huntingIntensity = NA)
huntingIntensity$huntingIntensity <- colSums(expDistances) 



################################################################################
######################### FORMAT SITECOV DATAFRAME #############################
################################################################################

#### edit site covariates ####
siteCovariate$HuntingIntensityRaw <- huntingIntensity$huntingIntensity
siteCovariate$HuntingIntensity <- c(scale(huntingIntensity$huntingIntensity))
siteCovariate$Station <- as.factor(paste0("ZAB", siteCovariate$Station))
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$CommunityRaw <- siteCovariate$Community/1000 # convert to km
siteCovariate$Community <- c(scale(siteCovariate$Community/1000))
siteCovariate$RiverRaw <- siteCovariate$River/1000 # convert to km
siteCovariate$River <- c(scale(siteCovariate$River/1000))
siteCovariate$EffortRaw <- siteCovariate$Effort
siteCovariate$Effort <- c(scale(siteCovariate$Effort))
siteCovariate$OnTrail <- as.factor(ifelse(siteCovariate$Trail.Distance == 0, 1, 0))
siteCovariate$Trail.DistanceRaw <- siteCovariate$Trail.Distance
siteCovariate$Trail.Distance <- c(scale(siteCovariate$Trail.Distance))
#siteCovariate$Station <- NULL
siteCovariate$RR <- NULL
siteCovariate$CR <- NULL


################################################################################
######################## FORMAT BUFFER COMPOSITION #############################
################################################################################
## How natural vs. non-natural are classified:
# (pulled from original DISES_BufferAnalysis Excel sheet by Michael W.)
# natural areas: forest formation, flooded forest, wetland, grassland, other non-forest natural area
# non-natural areas: farming, urban infrastructure, other non-vegetated area, mining

# make sure all the first columns match so we can use it as column names
all(all(ZABBuffer[,1]==SGEBuffer[,1]),
    all(SGEBuffer[,1]==SNABuffer[,1]),
    all(SNABuffer[,1]==SKPBuffer[,1]))
columns <- ZABBuffer[,1]

# transpose
ZABBufferT <- as.data.frame(t(ZABBuffer[,-1]))
colnames(ZABBufferT) <- columns

SGEBufferT <- as.data.frame(t(SGEBuffer[,-1]))
colnames(SGEBufferT) <- columns

SNABufferT <- as.data.frame(t(SNABuffer[,-1]))
colnames(SNABufferT) <- columns

SKPBufferT <- as.data.frame(t(SKPBuffer[,-1]))
colnames(SKPBufferT) <- columns

allBuffer <- rbind(ZABBufferT, SGEBufferT, SNABufferT, SKPBufferT)

# remove all rows with all NAs
allBuffer <- allBuffer[rowSums(is.na(allBuffer)) != ncol(allBuffer),]

# extract only necessary values
PercentNatural <- allBuffer[allBuffer$Year == "2020" & allBuffer$InOut == "Outside", 
                            c("Community", "PercentNaturalArea")]
names(PercentNatural)[names(PercentNatural) == 'Community'] <- 'CommunityName'

# merge the percent of natural area outside the community with all the stations
allRecords <- globalRecords[,c("CommunityName", "Station")]

################################################################################
################################## SAVE IT ##################################### 
################################################################################
write.csv(siteCovariate, file = 'siteCovs2018.csv')
write.csv(PercentNatural, file = '../../Global/Data/PercentNaturalOutside.csv')
write.csv(allBuffer, file = "../../Global/Data/allBuffer.csv")
save(siteCovariate, file = 'R Objects/siteCovs2018.RData')





























################################################################################
###################################### OLD ##################################### 
################################################################################

################ ASSIGNING EACH HUNTING OCCASION TO THE NEAREST CAMERA
# find the closest camera to each hunting occasion
closestCamera <- data.frame(huntingX = huntingCoords$x,
                            huntingY = huntingCoords$y,
                            distToNearestCam = NA,
                            station = NA)
for (i in 1:nrow(distMat)) {
  distances <- distMat[i,] # distance from hunt to all cameras
  closestCamera$station[i] <- stations$Station[which(distances == min(distances))]
  closestCamera$distToNearestCam[i] <- min(distances)
}

# only keep nearby hunts
threshold <- 2000 # in meters

includedHunts <- subset(closestCamera, distToNearestCam <= threshold)
excludedHunts <- subset(closestCamera, distToNearestCam > threshold)

# map the included hunts
includedHunts <- st_as_sf(includedHunts, coords = c("huntingX", "huntingY"),  crs = 4326)
excludedHunts <- st_as_sf(excludedHunts, coords = c("huntingX", "huntingY"),  crs = 4326)
#mapview(includedHunts, col.regions = "green", alpha = 0.1) +
#mapview(excludedHunts, col.regions = "red", alpha = 0.1) +
#mapview(cameras, alpha = 1)

#### results:
# within 500 meters = 40 hunts
# within 1000 meters = 176 hunts
# within 1500 meters = 249 hunts
# within 2000 meters = 297 hunts
# within 2500 meters = 386 hunts

################################################################## 