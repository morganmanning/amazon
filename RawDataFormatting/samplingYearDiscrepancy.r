
library(terra)
library(sf)
library(dplyr)
library(exactextractr)
library(ggplot2)
library(lme4)

setwd("~/Documents/amazon/Global/Data")

## camera trap data
cameras <- read.csv("AllStationsFormatted.csv")
coords <- data.frame(id = cameras[,c("Station")], 
lon = cameras[,c("gps_x")], 
lat = cameras[,c("gps_y")])

# buffer
bufferKM <- 10
cameraCRS <- "+proj=longlat +datum=WGS84 +no_defs" # cameras coming from: crs = 4326

# convert
points <- st_as_sf(coords, coords = c("lon", "lat"), crs = cameraCRS)
points_proj <- st_transform(points, cameraCRS)
buffers <- st_buffer(points_proj, bufferKM * 1000)

# list of files
files_2018 <- c(
  "GLDAS_NOAH025_M.A201807.021.nc4",
  "GLDAS_NOAH025_M.A201808.021.nc4",
  "GLDAS_NOAH025_M.A201809.021.nc4",
  "GLDAS_NOAH025_M.A201810.021.nc4",
  "GLDAS_NOAH025_M.A201811.021.nc4"
)

files_2022 <- c(
  "GLDAS_NOAH025_M.A202207.021.nc4",
  "GLDAS_NOAH025_M.A202208.021.nc4",
  "GLDAS_NOAH025_M.A202209.021.nc4",
  "GLDAS_NOAH025_M.A202210.021.nc4",
  "GLDAS_NOAH025_M.A202211.021.nc4"
)

# function to facilitate
extract_variable <- function(files, varname, buffers) {
  vals <- lapply(files, function(f) {
    r <- rast(f, subds = varname)
    exact_extract(r, buffers, 'mean')
  })
  
  do.call(cbind, vals) %>% as.data.frame()
}

# extract temp
temp_2018 <- extract_variable(paste0("~/Downloads/", files_2018), "Tair_f_inst", buffers)
temp_2022 <- extract_variable(paste0("~/Downloads/", files_2022), "Tair_f_inst", buffers)
temp_2018 <- temp_2018 - 273.15 # to convert from K to C
temp_2022 <- temp_2022 - 273.15 

# extract precip
rain_2018 <- extract_variable(paste0("~/Downloads/", files_2018), "Rainf_f_tavg", buffers)
rain_2022 <- extract_variable(paste0("~/Downloads/", files_2022), "Rainf_f_tavg", buffers)
seconds_per_month <- 30 * 24 * 60 * 60  # approx
rain_2018 <- rain_2018 * seconds_per_month
rain_2022 <- rain_2022 * seconds_per_month

# summary df
summary_df <- data.frame(
  id = coords$id,
  temp_2018 = rowMeans(temp_2018, na.rm = TRUE),
  temp_2022 = rowMeans(temp_2022, na.rm = TRUE),
  rain_2018 = rowMeans(rain_2018, na.rm = TRUE),
  rain_2022 = rowMeans(rain_2022, na.rm = TRUE)
)

print(summary_df)

# test it
t_temp <- t.test(summary_df$temp_2018, summary_df$temp_2022, paired = TRUE)
t_rain <- t.test(summary_df$rain_2018, summary_df$rain_2022, paired = TRUE)

# for smaller n
w_temp <- wilcox.test(summary_df$temp_2018, summary_df$temp_2022, paired = TRUE)
w_rain <- wilcox.test(summary_df$rain_2018, summary_df$rain_2022, paired = TRUE)

# more
long_df <- data.frame(
  site = rep(coords$id, each = 5 * 2),
  year = rep(c(2018, 2022), each = 5 * nrow(coords)),
  month = rep(7:11, times = 2 * nrow(coords)),
  temp = c(as.vector(as.matrix(temp_2018)),
           as.vector(as.matrix(temp_2022))),
  rain = c(as.vector(as.matrix(rain_2018)),
           as.vector(as.matrix(rain_2022)))
)

# all results:
w_temp
w_rain
t_temp
t_rain
lmer(temp ~ year + (1 | site), data = long_df)
lmer(rain ~ year + (1 | site), data = long_df)



###################### NATURAL AREA NEXT #############################
## 2018
r17M2018 <- rast("~/Downloads/17M_20180101-20190101.tif")
r17M2018_p <- project(r17M2018, cameraCRS, method = "near")
r17N2018 <- rast("~/Downloads/17N_20180101-20190101.tif")
r17N2018_p <- project(r17N2018, cameraCRS, method = "near")
r18M2018 <- rast("~/Downloads/18M_20180101-20190101.tif")
r18M2018_p <- project(r18M2018, cameraCRS, method = "near")
r18N2018 <- rast("~/Downloads/18N_20180101-20190101.tif")
r18N2018_p <- project(r18N2018, cameraCRS, method = "near")

## 2022
r17M2022 <- rast("~/Downloads/17M_20220101-20230101.tif")
r17M2022_p <- project(r17M2022, cameraCRS, method = "near")
r17N2022 <- rast("~/Downloads/17N_20220101-20230101.tif")
r17N2022_p <- project(r17N2022, cameraCRS, method = "near")
r18M2022 <- rast("~/Downloads/18M_20220101-20230101.tif")
r18M2022_p <- project(r18M2022, cameraCRS, method = "near")
r18N2022 <- rast("~/Downloads/18N_20220101-20230101.tif")
r18N2022_p <- project(r18N2022, cameraCRS, method = "near")

# merge tiles per year
# 2018
lulc_2018 <- mosaic(r17M2018_p, r17N2018_p, r18M2018_p, r18N2018_p)
writeRaster(lulc_2018, "lulc_2018.tif", overwrite = TRUE)
#plot(lulc_2018)

# 2022
lulc_2022 <- mosaic(r17M2022_p, r17N2022_p, r18M2022_p, r18N2022_p)
writeRaster(lulc_2022, "lulc_2022.tif", overwrite = TRUE)
#plot(lulc_2022)

# aggregate
agg_2018 <- aggregate(lulc_2018, fact = 3, fun = "modal")
writeRaster(agg_2018, "agg_lulc_2018.tif", overwrite = TRUE)
agg_2022 <- aggregate(lulc_2022, fact = 3, fun = "modal")
writeRaster(agg_2022, "agg_lulc_2022.tif", overwrite = TRUE)

# buffer
sitesBuffered <- st_buffer(
  st_as_sf(cameras[,c("gps_x", "gps_y")],
           coords = c("gps_x", "gps_y"),
           crs = 4326),
  bufferKM * 1000
)

# extract function
sum_cover <- function(x){
  list(x %>%
         mutate(across('value', round)) %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)) %>%
         mutate(proportion = total_area / sum(total_area)))
}

# extract both years
# 2018
x2018 <- exact_extract(agg_2018, sitesBuffered,
                       coverage_area = TRUE,
                       summarize_df = TRUE,
                       fun = sum_cover)

names(x2018) <- cameras$Station
lulc2018_df <- bind_rows(x2018, .id = "Station")

# 2022
x2022 <- exact_extract(agg_2022, sitesBuffered,
                       coverage_area = TRUE,
                       summarize_df = TRUE,
                       fun = sum_cover)

names(x2022) <- cameras$Station
lulc2022_df <- bind_rows(x2022, .id = "Station")

# refine
natural_vals <- c(1, 2, 4)  # water, trees, flooded vegetation

# nat per site
calc_percent_natural <- function(df) {
  df$value <- as.integer(round(df$value))
  
  df %>%
    group_by(Station) %>%
    summarize(percentNatural = sum(proportion[value %in% natural_vals], na.rm = TRUE))
}

nat2018 <- calc_percent_natural(lulc2018_df)
nat2022 <- calc_percent_natural(lulc2022_df)

# merge and compute change
nat_compare <- merge(nat2018, nat2022, by = "Station",
                     suffixes = c("_2018", "_2022"))

nat_compare$change <- nat_compare$percentNatural_2022 - nat_compare$percentNatural_2018

print(nat_compare)

# paired t-test
t_nat <- t.test(nat_compare$percentNatural_2018,
                nat_compare$percentNatural_2022,
                paired = TRUE)

t_nat

# smaller sample size safe
wilcox.test(nat_compare$percentNatural_2018,
            nat_compare$percentNatural_2022,
            paired = TRUE)
