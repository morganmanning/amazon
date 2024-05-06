##### formatting .tif land use land cover .tif files
# PROCESS: downloaded LULC data at a 10m resolution for the entire area
# PURPOSE: cut down these .tif files because they're too big to put on GitHub right now
# Coordinate system of TIFs: WGS84; UTM; ESPG:3857

############ SET UP ############ 
# set working directory
setwd("/Users/morganmanning/Documents/amazon/Global")
setwd("~/Documents/amazon/Global")

# load packages
require(terra)
require(sf)
require(exactextractr)
require(dplyr)

# load in TIFF files: https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/download.html (2020)
N00_W80 <- rast("../../../Downloads/00N_080W.tif")
N00_W90 <- rast("../../../Downloads/00N_090W.tif")
N10_W80 <- rast("../../../Downloads/10N_080W.tif")


## camera trap data
cameras <- read.csv("Data/AllStationsFormatted.csv")

################ COMBINE AND CUT THE RASTERS ##################
# make the camera coordinates into a matrix and project them
cameraCRS <- "+proj=longlat +datum=WGS84 +no_defs" # cameras coming from: crs = 4326
camCoordMatrix <- as.matrix(cameras[,c("gps_x", "gps_y")])

# make sure matching projections
tifCRS <- crs(N00_W80, proj = TRUE)

# buffered extent
bufferDegrees <- 1.5
e <- as.vector(ext(camCoordMatrix)) 
e["xmin"] <- e["xmin"] - bufferDegrees
e["ymin"] <- e["ymin"] - bufferDegrees
e["xmax"] <- e["xmax"] + bufferDegrees
e["ymax"] <- e["ymax"] + bufferDegrees
eProjected <- ext(e)
plot(eProjected)

# crop them
N00_W80_cropped <- crop(N00_W80, eProjected)
# N00_W90_cropped <- crop(N00_W90, eProjected) # doesn't overlap
N10_W80_cropped <- crop(N10_W80, eProjected)

# merge them
raster2020 <- merge(N00_W80_cropped, N10_W80_cropped)
plot(N00_W80_cropped)
plot(N10_W80_cropped)
plot(raster2020)

# save it
writeRaster(raster2020, "Data/cropped2020raster.tif", overwrite = TRUE)
raster2020 <- rast("Data/cropped2020raster.tif")

################ EXTRACT LULC FROM BUFFERED SITES ##################

# buffer the points 
bufferKM <- 25
bufferedPoints <- buffer(vect(camCoordMatrix, type = "points", crs = cameraCRS), width = bufferKM*1000)
plot(bufferedPoints)

# extract proportion of each land cover per each buffered site
sites <- st_as_sf(cameras[,c("gps_x", "gps_y")], coords = c("gps_x", "gps_y"), crs = cameraCRS)
sitesBuffered <- st_buffer(sites, bufferKM*1000)

sum_cover <- function(x){
  list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)) %>%
         mutate(proportion = total_area/sum(total_area)))
  
}

# extract the area of each raster cell covered by the plot and summarize
x <- exact_extract(raster2020, sitesBuffered, coverage_area = TRUE, 
                   summarize_df = TRUE, fun = sum_cover)

# add plot names to the elements of the output list
names(x) <- cameras$Station

# merge the list elements into a df
LULCperSite <- bind_rows(x, .id = "Station")





















##################### BELOW WOULD NOT WORK
# SOURCE: https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d; 
# DOWNLOAD: https://livingatlas.arcgis.com/landcoverexplorer/#mapCenter=39.18600%2C9.04200%2C10&mode=step&timeExtent=2017%2C2022&year=2020&downloadMode=true 

# some of the .tifs had the same extents, making it impossible to crop
# when I plotted the the extent box or points on top of the rasters that cropped, neither were visible

#################################
# read in all .tif files 
## 2018
r17M2018 <- rast("../../../Downloads/17M_20180101-20190101.tif")
r17N2018 <- rast("../../../Downloads/17N_20180101-20190101.tif")
r18M2018 <- rast("../../../Downloads/18M_20180101-20190101.tif")
r18N2018 <- rast("../../../Downloads/18N_20180101-20190101.tif")

## 2022
r17M2022 <- rast("../../../Downloads/17M_20220101-20230101.tif")
r17N2022 <- rast("../../../Downloads/17N_20220101-20230101.tif")
r18M2022 <- rast("../../../Downloads/18M_20220101-20230101.tif")
r18N2022 <- rast("../../../Downloads/18N_20220101-20230101.tif")

require(raster)
## 2020
M18 <- rast("../../../Downloads/18M_20200101-20210101.tif")
M17 <- rast("../../../Downloads/17M_20200101-20210101.tif")
N18 <- rast("../../../Downloads/18N_20200101-20210101.tif")
N17 <- rast("../../../Downloads/17N_20200101-20210101.tif")

##### the issue
par(mfrow = c(1,2))
plot(M18) ; plot(M17)
ext(M18) == ext(M17)

## tried:
# - using 'raster' package instead of 'terra'
# - using different rasters (it just had land cover, not land use)














## camera trap data
cameras <- read.csv("Data/AllStationsFormatted.csv")

################ COMBINE AND CUT THE RASTERS ##################
# make the camera coordinates into a matrix and project them
cameraCRS <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" # cameras coming from: crs = 4326
camerasSV <- vect(as.matrix(cameras[,c("gps_x", "gps_y")]), crs = cameraCRS)

# buffered extent
bufferDegrees <- 1.5
e <- as.vector(ext(camerasSV)) 
e["xmin"] <- e["xmin"] - bufferDegrees
e["ymin"] <- e["ymin"] - bufferDegrees
e["xmax"] <- e["xmax"] + bufferDegrees
e["ymax"] <- e["ymax"] + bufferDegrees
eBuffered <- ext(e)
ext(camerasSV); eBuffered

# project camera extent into UTM to crop
eUTM <- project(ext(eBuffered), from = cameraCRS, to = crs(M18, proj = TRUE))
eUTM

# crop
M18_cropped <- crop(M18, eUTM)
M17_cropped <- crop(M17, eUTM)
# N18_cropped <- crop(N18, eUTM)
# N17_cropped <- crop(N17, eUTM)


# project tifs
M18_latlon <- project(M18_cropped, y = "epsg:4326")
M17_latlon <- project(M17_cropped, y = "epsg:4326")
N18_latlon <- project(N18_cropped, y = "epsg:4326")
N17_latlon <- project(N17_cropped, y = "epsg:4326")



merged <- merge(M18, M17)
merged <- merge(merged, N18)
merged <- merge(merged, N17)




camCoordMatrix <- as.matrix(cameras[,c("gps_x", "gps_y")])
projectedPoints <- project(camCoordMatrix, from = cameraCRS, to = tif18CRS)

# buffered extent
bufferKM <- 200
e <- as.vector(ext(projectedPoints)) 
e["xmin"] <- e["xmin"] - (1000*bufferKM) 
e["ymin"] <- e["ymin"] - (1000*bufferKM) 
e["xmax"] <- e["xmax"] + (1000*bufferKM) 
e["ymax"] <- e["ymax"] + (1000*bufferKM) 
eProjected <- ext(e)

# crop them
# 2018
cr17M2018 <- crop(r17M2018, eProjected)
# cr17N2018 <- crop(r17N2018, eProjected) # Error: [crop] extents do not overlap
cr18M2018 <- crop(r18M2018, eProjected)
# cr18N2018 <- crop(r18N2018, eProjected) # Error: [crop] extents do not overlap

# 2022
cr17M2022 <- crop(r17M2022, eProjected)
# cr17N2022 <- crop(r17N2022, eProjected) # Error: [crop] extents do not overlap
cr18M2022 <- crop(r18M2022, eProjected)
# cr18N2022 <- crop(r18N2022, eProjected) # Error: [crop] extents do not overlap

# combine them by year
raster2018 <- mosaic(cr17M2018, cr18M2018)
raster2022 <- mosaic(cr17M2022, cr18M2022)


### something is wrong because these don't look the same
plot(raster2018)
plot(projectedPoints)

# save them
#writeRaster(raster2018, "Data/cropped2018raster.tif", overwrite = TRUE)
#writeRaster(raster2022, "Data/cropped2022raster.tif", overwrite = TRUE)



