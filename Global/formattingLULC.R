##### formatting .tif land use land cover .tif files
# SOURCE: https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d
# PROCESS: downloaded LULC data at a 10m resolution for the entire area
# PURPOSE: cut down these .tif files because they're too big to put on GitHub right now
# Coordinate system of TIFs: WGS84; UTM; ESPG:3857

############ SET UP ############ 
# set working directory
setwd("/Users/morganmanning/Documents/amazon/Global")

# load packages
require(terra)

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

## camera trap data
cameras <- read.csv("Data/AllStationsFormatted.csv")

################ COMBINE AND CUT THE RASTERS ##################
# make the camera coordinates into a matrix and project them
cameraCRS <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" # cameras coming from: crs = 4326
tifCRS <- crs(r17M2018, proj = TRUE)
camCoordMatrix <- as.matrix(cameras[,c("gps_x", "gps_y")])
projectedPoints <- project(camCoordMatrix, from = cameraCRS, to = tifCRS)

# buffered extent
e <- as.vector(ext(projectedPoints)) 
e["xmin"] <- e["xmin"] - (1000*100) # buffer by 100km
e["ymin"] <- e["ymin"] - (1000*100) # buffer by 100km
e["xmax"] <- e["xmax"] + (1000*100) # buffer by 100km
e["ymax"] <- e["ymax"] + (1000*100) # buffer by 100km
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
raster2018 <- merge(cr17M2018, cr18M2018)
raster2022 <- merge(cr17M2022, cr18M2022)


### something is wrong because these don't look the same
plot(raster2018)
points(projectedPoints)



# save them
writeRaster(raster2018, "Data/cropped2018raster.tif", overwrite = TRUE)
writeRaster(raster2022, "Data/cropped2022raster.tif", overwrite = TRUE)



