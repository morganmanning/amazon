##### formatting .tif land use land cover .tif files
# SOURCE: https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d
# PROCESS: downloaded LULC data at a 10m resolution for the entire area
# PURPOSE: cut down these .tif files because they're too big to put on GitHub right now
# Coordinate system: WGS84; UTM; ESPG:3857

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
# combine them by year
raster2018 <- mosaic(r17M2018, r17N2018, r18M2018, r18N2018)
raster2022 <- mosaic(r17M2022, r17N2022, r18M2022, r18N2022)

# make a bounding box
maxX <- max(cameras$gps_x) + 1 # longitude
minX <- min(cameras$gps_x) - 1 # longitude
maxY <- max(cameras$gps_y) + 1 # latitude
minY <- min(cameras$gps_y) - 1 # latitude
e <- ext(minX, maxX, minY, maxY) # xmin, xmax, ymin, ymax

# crop each raster with the bounding box
cropped2018 <- crop(raster2018, e)
cropped2022 <- crop(raster2022, e)

# save them
writeRaster(cropped2018, "Data/cropped2018raster.tif", overwrite = TRUE)
writeRaster(cropped2022, "Data/cropped2022raster.tif", overwrite = TRUE)



