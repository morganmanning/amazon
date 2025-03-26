
# Goal: collect covariates for each community into one data frame for multidimensional scaling

setwd("~/Documents/amazon/Global/Data")

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d
# source for rainfall, temperature, humidity, root moisture: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary

################################################################################
################################################################################
################################################################################

# load packages
require(terra)
require(sf)
require(exactextractr)
require(dplyr)

# load necessary data
covariates <- read.csv("AllCommunityCovariates.csv")
stations <- read.csv("AllStationsFormatted.csv")
cameraCoords <- stations[, c("gps_x", "gps_y")]
Data <- read.csv("AllIndependentRecordsFormatted.csv")


# average the number of animal detections per day
nDetections <- Data %>%
    group_by(CommunityName) %>%
    mutate(nDetectionsTotal = n(), nDays = max(as.Date(DateTimeOriginal)) - min(as.Date(DateTimeOriginal))) %>%
    mutate(nDetectionsPerDay = as.numeric(nDetectionsTotal) / as.numeric(nDays))

nDetections <- nDetections[,c("CommunityName", "nDetectionsPerDay")]
nDetections <- nDetections[!duplicated(nDetections), ]

# rename "CommunityName" to "Community" for merging
colnames(nDetections) <- c("Community", "nDetectionsPerDay")

# add to covariates dataset
covariates <- merge(covariates, nDetections, by = "Community")



################################################################################
################################################################################
######################## ADD ELEVATION TO COVARIATES ###########################
################################################################################
################################################################################

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d

# load elevation data (only stored locally on work computer)
dem0 <- rast("Community-level Covariates/DEM/sa_dem_0.tif")
dem1 <- rast("Community-level Covariates/DEM/sa_dem_1.tif")
dem2 <- rast("Community-level Covariates/DEM/sa_dem_2.tif")
dem3 <- rast("Community-level Covariates/DEM/sa_dem_3.tif")
dem4 <- rast("Community-level Covariates/DEM/sa_dem_4.tif")

# merge rasters
all.equal(crs(dem0), crs(dem1), crs(dem2), crs(dem3), crs(dem4)) # make sure all in same crs
rasterCollection <- sprc(dem0, dem1, dem2, dem3, dem4)
mergedRaster <- merge(rasterCollection)

# extract values at each camera trap
stations$Elevation <- extract(mergedRaster, cameraCoords, xy = FALSE, ID = FALSE)$sa_dem_0
justElev <- stations[, c("Station", "Elevation")]
str(justElev)

# add elevation to covariates
covariates <- merge(justElev, covariates, by = "Station")



################################################################################
################################################################################
##################### ADD DISTANCE TO ANOTHER COMMUNITY ########################
################################################################################
################################################################################

# load community points
secoya <- read_sf("Community Points/Secoya.kml")
remolino <- secoya[secoya$Name == "Remolino", ]
sanPablo <- secoya[secoya$Name == "San Pablo", ]
sinangoe <- read_sf("Community Points/Sinangoe.kml")
siona <- read_sf("Community Points/Siona.kml")
zabalo <- read_sf("Community Points/Zabalo.kml")

# calculate distances between points
communityPoints <- data.frame(
    Community = c("Remolino", "San Pablo", "Sinangoe", "Siona", "Zabalo"),
    X = c(st_coordinates(remolino)[, "X"], st_coordinates(sanPablo)[, "X"], st_coordinates(sinangoe)[, "X"], colMeans(st_coordinates(siona))[c("X")], st_coordinates(zabalo)[, "X"]),
    Y = c(st_coordinates(remolino)[, "Y"], st_coordinates(sanPablo)[, "Y"], st_coordinates(sinangoe)[, "Y"], colMeans(st_coordinates(siona))[c("Y")], st_coordinates(zabalo)[, "Y"]),
    distanceToNearestCommunity = NA,
    TerritoryArea = c(212619183, 56178602, 307786455, 121874135, 1374033979) # gotten from putting territory shapefiles in qGIS and calculatin area
)
distances <- distance(communityPoints[, c("X", "Y")], communityPoints[, c("X", "Y")], lonlat = TRUE)
colnames(distances) <- communityPoints$Community
rownames(distances) <- communityPoints$Community
for (i in 1:nrow(communityPoints)) {
    # exclude 0 because that's the distance to itself
    dists <- distances[i, ]
    dists <- dists[dists != 0]
    communityPoints$distanceToNearestCommunity[i] <- min(dists)
}



################################################################################
################################################################################
##### ADD RAINFALL, TEMPERATURE, HUMIDITY, AND ROOT MOISTURE TO COVARIATES #####
################################################################################
################################################################################
# methodology for how I got these CSVs are in meeting_notes.md

# Load in the data
rainfall <- read.csv("Community-level Covariates/rainfall.csv")
humidity <- read.csv("Community-level Covariates/humidity.csv")
airTemp <- read.csv("Community-level Covariates/airTemp.csv")
rootMoisture <- read.csv("Community-level Covariates/rootMoisture.csv")

# Combine into one data frame
NASAcovariates <- rbind(rainfall, humidity, airTemp, rootMoisture)

colnames(NASAcovariates) <- gsub("Name", "Community", colnames(NASAcovariates))

# Change San Pablo/Siona North/South to San Pablo
NASAcovariates$Community <- gsub("San Pablo North", "San Pablo", NASAcovariates$Community)
NASAcovariates$Community <- gsub("San Pablo South", "San Pablo", NASAcovariates$Community)
NASAcovariates$Community <- gsub("Siona North", "Siona", NASAcovariates$Community)
NASAcovariates$Community <- gsub("Siona South", "Siona", NASAcovariates$Community)

# Average x mean and sd by community
NASAcovariates <- NASAcovariates %>%
    select(Community, Covariate, X_mean, X_stdev) %>%
    group_by(Community, Covariate) %>%
    summarize(
        Mean = mean(X_mean),
        SD = mean(X_stdev)
    )

# put it into the big data frame
communityPoints$rainfall <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "rainfall", ], by = "Community", all.x = TRUE)$Mean 
communityPoints$rainfallSD <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "rainfall", ], by = "Community", all.x = TRUE)$SD
communityPoints$humidity <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "humidity", ], by = "Community", all.x = TRUE)$Mean
communityPoints$humiditySD <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "humidity", ], by = "Community", all.x = TRUE)$SD
communityPoints$airTemp <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "airTemp", ], by = "Community", all.x = TRUE)$Mean
communityPoints$airTempSD <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "airTemp", ], by = "Community", all.x = TRUE)$SD
communityPoints$rootMoisture <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "rootMoist", ], by = "Community", all.x = TRUE)$Mean
communityPoints$rootMoistureSD <- merge(communityPoints, NASAcovariates[NASAcovariates$Covariate == "rootMoist", ], by = "Community", all.x = TRUE)$SD



################################################################################
################################################################################
##################### ADD SPECIES DIVERSITY AND RICHNESS #######################
################################################################################
################################################################################

# load species data
diversity <- read.csv("CommunityDiversityAbundance.csv")
diversity$Community <- gsub("Zábalo", "Zabalo", diversity$Community)
diversity$X <- NULL

# merge community points with diversity data
communityCovariates <- merge(communityPoints, diversity, by = "Community")
communityCovariates$PopulationSize <- c(NA, NA, 104, 116, 87) # taken from AAG paper table #1
communityCovariates$DaysHuntingPerMonthDry <- c(NA, NA, 2.37, 1.36, 4.62)
communityCovariates$DaysHuntingPerMonthWet <- c(NA, NA, 0.84, 0.43, 1.98)
communityCovariates$DaysFishingPerMonthDry <- c(NA, NA, 2.27, 2.38, 5.30)
communityCovariates$DaysFishingPerMonthWet <- c(NA, NA, 1.22, 0.46, 1.10)
communityCovariates$PercentPopWhoHunt <- c(NA, NA, 67, 74, 49)
communityCovariates$PercentPopWhoFish <- c(NA, NA, 80, 95, 98)

meanPhysical <- covariates %>%
    group_by(Community) %>%
    summarize(
        # MeanRainfall = mean(Rainfall), # Rainfall
        MeanTemperature = mean(Temperature),
        MeanDistToWater = mean(DistToWater),
        MeanElevation = mean(Elevation)
    )
communityCovariates <- merge(communityCovariates, meanPhysical, by = "Community")

# save as a csv
write.csv(communityCovariates, "CommunityLevelCovariates.csv", row.names = FALSE)
