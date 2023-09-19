
######################## SETTING UP SITE COVARIATES ##########################
# set my working directory
# setwd("~/Documents/amazon/Data")

# load required packages
require(ggplot2)
require(dplyr)
require(sf)
require(mapview)
require(geosphere)

# load files
siteCovariate <- read.csv("siteCovs2018Raw.csv")
stations <- read.csv("Stations2018Raw.csv")
stations <- stations %>% 
  distinct(x, y, .keep_all = TRUE)
trapRecords <- read.csv("RecordTable2018.csv") # trapped species records
hunting <- read.csv("HuntingData2018.csv") # harvest species, location, and date




################################################################################
######################## FORMAT HUNTING INTENSITY ##############################
################################################################################

# formatting the data as spatial objects
# latitude = y, longitude = x
cameras <- st_as_sf(stations, coords = c("x", "y"),  crs = 4326) # EPSG:24878 = UTM zone 18S
harvests <- st_as_sf(hunting, coords = c("x", "y"),  crs = 4326)
mapview(harvests, col.regions = "red", alpha = 0.1) +
  mapview(cameras, alpha = 1)

# just the coordinates
cameraCoords <- stations[,c('x', 'y')]
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
siteCovariate$HuntingIntensity <- scale(huntingIntensity$huntingIntensity)
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$Community <- scale(siteCovariate$Community/1000)
siteCovariate$River <- scale(siteCovariate$River/1000)
siteCovariate$Effort <- scale(siteCovariate$Effort)
siteCovariate$OnTrail <- as.factor(ifelse(siteCovariate$Trail.Distance == 0, 1, 0))
siteCovariate$Trail.Distance <- scale(siteCovariate$Trail.Distance)
siteCovariate$Station <- NULL
siteCovariate$RR <- NULL
siteCovariate$CR <- NULL





################################################################################
################################## SAVE IT ##################################### 
################################################################################
write.csv(siteCovariate, file = 'siteCovs2018.csv')
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