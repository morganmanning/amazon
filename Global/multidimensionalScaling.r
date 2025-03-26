# Goal: collect covariates for each community into one data frame for multidimensional scaling

setwd("~/Documents/amazon/Global/Data")

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d

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


################################################################################
################################################################################
######################## ADD ELEVATION TO COVARIATES ###########################
################################################################################
################################################################################

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d

# load elevation data
dem0 <- rast("DEM/sa_dem_0.tif")
dem1 <- rast("DEM/sa_dem_1.tif")
dem2 <- rast("DEM/sa_dem_2.tif")
dem3 <- rast("DEM/sa_dem_3.tif")
dem4 <- rast("DEM/sa_dem_4.tif")

# merge rasters
all.equal(crs(dem0), crs(dem1), crs(dem2), crs(dem3), crs(dem4)) # make sure all in same crs
rasterCollection <- sprc(dem0, dem1, dem2, dem3, dem4)
mergedRaster <- merge(rasterCollection)

# extract values at each camera trap
stations$Elevation <- extract(mergedRaster, cameraCoords, xy = FALSE, ID = FALSE)$sa_dem_0
justElev <- stations[,c("Station", "Elevation")]
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
remolino <- secoya[secoya$Name == "Remolino",]
sanPablo <- secoya[secoya$Name == "San Pablo",]
sinangoe <- read_sf("Community Points/Sinangoe.kml")
siona <- read_sf("Community Points/Siona.kml")
zabalo <- read_sf("Community Points/Zabalo.kml")

# calculate distances between points
communityPoints <- data.frame(Community = c("Remolino", "San Pablo", "Sinangoe", "Siona", "Zabalo"),
                                X = c(st_coordinates(remolino)[,"X"], st_coordinates(sanPablo)[,"X"], st_coordinates(sinangoe)[,"X"], colMeans(st_coordinates(siona))[c("X")], st_coordinates(zabalo)[,"X"]),
                                Y = c(st_coordinates(remolino)[,"Y"], st_coordinates(sanPablo)[,"Y"], st_coordinates(sinangoe)[,"Y"], colMeans(st_coordinates(siona))[c("Y")], st_coordinates(zabalo)[,"Y"]),
                                distanceToNearestCommunity = NA, 
                                TerritoryArea = c(212619183, 56178602, 307786455, 121874135, 1374033979) # gotten from putting territory shapefiles in qGIS and calculatin area
                                )
distances <- distance(communityPoints[,c("X", "Y")], communityPoints[,c("X", "Y")], lonlat = TRUE)
colnames(distances) <- communityPoints$Community
rownames(distances) <- communityPoints$Community
for (i in 1:nrow(communityPoints)) {
    # exclude 0 because that's the distance to itself
    dists <- distances[i,]
    dists <- dists[dists != 0]
    communityPoints$distanceToNearestCommunity[i] <- min(dists)
}



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
communityCovariates$DaysFishingPerMonthWet <- c(NA, NA, 1.22, 0.46 , 1.10)
communityCovariates$PercentPopWhoHunt <- c(NA, NA, 67, 74, 49)
communityCovariates$PercentPopWhoFish <- c(NA, NA, 80, 95, 98)

meanPhysical <- covariates %>%
    group_by(Community) %>%
    summarize(
        MeanRainfall = mean(Rainfall), 
        MeanTemperature = mean(Temperature),
        MeanDistToWater = mean(DistToWater),
        MeanElevation = mean(Elevation)
    )
communityCovariates <- merge(communityCovariates, meanPhysical, by = "Community")
