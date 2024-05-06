##### formatting .tif land use land cover .tif files
# PROCESS: downloaded LULC data at a 10m resolution for the entire area
# PURPOSE: cut down these .tif files because they're too big to put on GitHub right now
# Coordinate system of TIFs: WGS84; UTM; ESPG:3857

############ SET UP ############ 
# set working directory
setwd("/Users/morganmanning/Documents/amazon/Global")

# load packages
require(terra)

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











##################### BELOW WOULD NOT WORK
# SOURCE: https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d; 

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

## camera trap data
cameras <- read.csv("Data/AllStationsFormatted.csv")

################ COMBINE AND CUT THE RASTERS ##################
# make the camera coordinates into a matrix and project them
cameraCRS <- "+proj=longlat +datum=WGS84 +no_defs +type=crs" # cameras coming from: crs = 4326
tif17CRS <- crs(r17M2018, proj = TRUE)
tif18CRS <- crs(r18M2018, proj = TRUE)

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
writeRaster(raster2018, "Data/cropped2018raster.tif", overwrite = TRUE)
writeRaster(raster2022, "Data/cropped2022raster.tif", overwrite = TRUE)



