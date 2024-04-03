
###############################################################################
###############################################################################
############### PLOTTING HUNTING RECORDS VS. CAMERA LOCATIONS #################
###############################################################################
###############################################################################



################################ GET SET UP ###################################
# set my working directory
# setwd("~/Documents/amazon/Data")

# load required packages
require(ggplot2)
require(dplyr)
require(sf)
require(mapview)
require(geosphere)

# load in data
stations <- read.csv("Stations2018.csv") # camera locations and dates
stations <- stations %>% 
  distinct(x, y, .keep_all = TRUE)
trapRecords <- read.csv("RecordTable2018.csv") # trapped species records
hunting <- read.csv("HuntingData2018.csv") # harvest species, location, and date




####################### PLOT TRAP AND HUNT LOCATIONS ###########################
# latitude = y, longitude = x
cameras <- st_as_sf(stations, coords = c("x", "y"),  crs = 4326) # EPSG:24878 = UTM zone 18S
harvests <- st_as_sf(hunting, coords = c("x", "y"),  crs = 4326)
mapview(harvests, col.regions = "red", alpha = 0.1) +
  mapview(cameras, alpha = 1)

cameraCoords <- stations[,c('x', 'y')]
huntingCoords <- hunting[,c('x', 'y')]

# distance from each trap to each hunting occasion
distMat <- distm(as.matrix(huntingCoords), as.matrix(cameraCoords), fun = distGeo)


####################### ASSIGN EACH HUNT TO A CAMERA ###########################
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


########################## EDIT CAMERA ASSIGNMENTS #############################
# several of the hunts are far away from any camera, so I'm going to remove far hunts

threshold <- 2000 # in meters

includedHunts <- subset(closestCamera, distToNearestCam <= threshold)
nrow(includedHunts)
excludedHunts <- subset(closestCamera, distToNearestCam > threshold)

# map the included hunts
includedHunts <- st_as_sf(includedHunts, coords = c("huntingX", "huntingY"),  crs = 4326)
excludedHunts <- st_as_sf(excludedHunts, coords = c("huntingX", "huntingY"),  crs = 4326)
mapview(includedHunts, col.regions = "green", alpha = 0.1) +
  mapview(excludedHunts, col.regions = "red", alpha = 0.1) +
  mapview(cameras, alpha = 1)

#### results:
# within 500 meters = 40 hunts
# within 1000 meters = 176 hunts
# within 1500 meters = 249 hunts
# within 2000 meters = 297 hunts
# within 2500 meters = 386 hunts


#################### WEIGHT HUNTS FOR EACH CAMERA TRAP #########################
### goal: value of hunting intensity for each camera trap

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

#### save it
save(huntingIntensity, file = "R Objects/huntingIntensity.RData")







# include all the hunting occasions, but weight the closer sites higher than the further sites 
# HuntingIntensity(site_i) = SIGMA(e^-(distance between site_i and hunting occasion_h)) 
# from h = 777 (above sigma) to h = 1 (below sigma)
# output will be a list of hunting intensity for each camera trap



# Jessica Kahler does wildlife crime
# CCJ 5934

# back to one-pager
# poaching cybercrime in the Amazon
# questions on the second one
# put one pagers in GoogleDocs 








###### plotting theoretical plots #######
# distance to community versus occupancy probability

# occupancy probability 
occProb <- c(0.01, 0.038, 0.05, 0.1, 0.25, 0.5, 0.8, 0.95, 0.99)
distanceToComm <- c(seq(0, 50, length.out = length(occProb)))

occProb <- c(0.99, 0.95, 0.8, 0.5, 0.25, 0.1, 0.05, 0.038, 0.01)
huntingFreq <- seq(0, 10, length.out = length(occProb))


####
plot(occProb ~ distanceToComm, type = 'b')

df <- data.frame(occProb, distanceToComm)
ggplot(df, aes(distanceToComm, occProb)) +
  geom_smooth(color = "black", se = FALSE) + 
  geom_ribbon(aes(ymin = occProb - 0.1, ymax = occProb + 0.1), alpha = 0.3) +
  #xlim(0, max(distanceToComm)) +
  #ylim(0, 1) +
  xlab("Distance to the community (km)") +
  ylab("Occupancy probability")

df <- data.frame(occProb, huntingFreq)
ggplot(df, aes(huntingFreq, occProb)) +
  geom_smooth(color = "black", se = FALSE) + 
  geom_ribbon(aes(ymin = occProb - 0.1, ymax = occProb + 0.1), alpha = 0.3) +
  #xlim(0, max(distanceToComm)) +
  #ylim(0, 1) +
  xlab("Hunting frequency (hunts/week)") +
  ylab("Occupancy probability")






