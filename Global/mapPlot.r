
# set working directory
setwd("~/Documents/amazon")

# load packages
require(dplyr)
require(ggplot2)
require(reshape2)
require(rphylopic)
require(ggpubr)
require(sf)
require(spData)
require(ggmagnify)
require(terra)
require(mapdata)
require(knitr)
require(lubridate)
require(tidyverse)
require(rnaturalearth)
require(ggrepel)
require(cowplot)
require(ggspatial) # for scale bar and north arrow

################################################################################
########################### PLOT MAP OF STATIONS ###############################
################################################################################

# load data
data("world")
stations <- read.csv("Global/Data/AllStationsFormatted.csv")
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")

# format data
stations <- stations %>%
  select(c(Station, gps_x, gps_y, CommunityName)) %>%
  distinct()
stations$CommunityName <- gsub("Zabalo", "Zábalo", x = stations$CommunityName)

# order the communities by percent of natural cover
orderedCommunities <- c(siteCovariate %>%
    group_by(Community) %>%
    summarize(perc = mean(PercentNaturalScaled)) %>%
    arrange(desc(perc)) %>%
    select(Community))$Community
orderedCommunities <- gsub("Zabalo", "Zábalo", x = orderedCommunities)

# order the communities by percent of natural cover
stations <- stations %>% dplyr::filter(CommunityName %in% orderedCommunities)
stations$Community <- factor(stations$CommunityName, levels=orderedCommunities)
stations$CommunityName <- NULL

# only highlight Ecuador
SA <- c("ecuador", "bolivia", "brazil", "chile", "colombia", "argentina", "guyana", "paraguay", "peru", "suriname", "uruguay", "venezuela")
mapColors <- rep("white", length(SA))
mapColors[2] <- "lightyellow"
colors <- c(
    "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
)

# higher-resolution South America
south_america_sf <- ne_countries(
    scale = "large",
    continent = "south america",
    returnclass = "sf"
)
south_america_sf$Ecu <- ifelse(south_america_sf$name_en == "Ecuador", "A", "B")

# import .kml community points
remSP <- sf::st_read("../One Drive Copy and Covariates/Secoya.kml")
remSP <- remSP %>%
    # dplyr::filter(Name == "Remolino" | Name == "San Pablo") %>%
    dplyr::mutate(Community = Name)

siona <- sf::st_read("../One Drive Copy and Covariates/Siona.kml")
siona <- siona %>%
    dplyr::mutate(Community = "Siona")

sinangoe <- sf::st_read("../One Drive Copy and Covariates/Sinangoe.kml")
sinangoe <- sinangoe %>%
    dplyr::mutate(Community = Name)

zabalo <- sf::st_read("../One Drive Copy and Covariates/Zabalo.kml")
zabalo <- zabalo %>%
    dplyr::mutate(Community = Name) %>%
    dplyr::mutate(Community = ifelse(Community == "Cofan Zabalo", "Zábalo", Community))

# all communities together
communities_sf <- rbind(remSP, siona, sinangoe, zabalo)

# load in territory .shp files
sanPabloTerritory <- st_read("../One Drive Copy and Covariates/Territories/SanPablo/SanPablo.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "San Pablo")
remolinoTerritory <- st_read("../One Drive Copy and Covariates/Territories/Remolino/Remolino.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Remolino")
siekopaiTerritory <- st_read("../One Drive Copy and Covariates/Territories/Siekopai/Siekopai.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry))
sionaTerritory <- st_read("../One Drive Copy and Covariates/Territories/Siona/Siona.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Siona")
sinangoeTerritory <- st_read("../One Drive Copy and Covariates/Territories/Sinangoe/Sinangoe.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Sinangoe")
zabaloTerritory <- st_read("../One Drive Copy and Covariates/Territories/Zabalo/Zabalo.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Zábalo")
# combine territories
territories_sf <- rbind(
    sanPabloTerritory, remolinoTerritory,
    sionaTerritory,
    sinangoeTerritory, zabaloTerritory
)

# set bounding box to magnify
coords <- as.matrix(stations[, c("gps_x", "gps_y")])
e <- as.vector(ext(coords))
e["xmin"] <- e["xmin"] - 0.1
e["ymin"] <- e["ymin"] - 0.5
e["xmax"] <- e["xmax"] + 0.1
e["ymax"] <- e["ymax"] + 0.5
# OG_inlay <- c(xmin = -150, xmax = -85, ymin = -45, ymax = 3)
new_inlay <- c(xmin = -190, xmax = -85, ymin = -40, ymax = 40)

# plot it
ecuadorMap <- south_america_sf %>%
    dplyr::filter(continent == "South America") %>%
    ggplot() +
    geom_sf(aes(fill = Ecu), lwd = 0.5) +
    geom_sf(data = territories_sf, aes(fill = Community), alpha = 0.4, lwd = 0.5) +
    geom_point(data = stations, aes(gps_x, gps_y, fill = Community), pch = 21, size = 1) +
    geom_sf(data = communities_sf, aes(fill = Community), size = 4, pch = 24) +
    geom_text_repel(
        data = communities_sf,
        aes(label = Community, geometry = geometry),
        box.padding = 0.75,
        max.overlaps = Inf,
        family = "Times",
        stat = "sf_coordinates",
        seed = 123 # for reproducible results
    ) +
    scale_fill_manual(
        name = "Community",
        values = c(colors, "A" = "lightyellow", "B" = "white"),
        breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    ) +
    coord_sf(
        default_crs = sf::st_crs(4326),
        xlim = c(new_inlay["xmin"], -37),
        ylim = c(-55, new_inlay["ymax"])
    ) +
    geom_magnify(
        from = e,
        to = new_inlay,
        shadow = TRUE
    ) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black")
    )
ecuadorMap




