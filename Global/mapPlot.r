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
require(tidyterra) # for geom_spatraster
require(ggnewscale)


################################################################################
########################### PLOT MAP OF STATIONS ###############################
################################################################################
# load data
data("world")
stations <- read.csv("Global/Data/AllStationsFormatted.csv")
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
naturalArea <- rast("Global/Data/Community-level Covariates/DEM/naturalArea.tif")
naturalArea <- as.factor(naturalArea)

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
stations$Community <- factor(stations$CommunityName, levels = orderedCommunities)
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
remSP <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Secoya.kml")
remSP <- remSP %>%
    # dplyr::filter(Name == "Remolino" | Name == "San Pablo") %>%
    dplyr::mutate(Community = ifelse(Name == "San Pablo" | Name == "Waiya" | Name == "Bellavista", "San Pablo", "Remolino"))
siona <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Siona.kml")
siona <- siona %>%
    dplyr::mutate(Community = "Siona")
sinangoe <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Sinangoe.kml")
sinangoe <- sinangoe %>%
    dplyr::mutate(Community = "Sinangoe")
zabalo <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Zabalo.kml")
zabalo <- zabalo %>%
    dplyr::mutate(Community = "Zábalo") %>%
    dplyr::mutate(Name = ifelse(Name == "Cofan Zabalo", "Zábalo", Name))

# all communities together
communities_sf <- rbind(remSP, siona, sinangoe, zabalo)

# load in territory .shp files
sanPabloTerritory <- st_read("Global/Data/Community-level Covariates/Territories/SanPablo/SanPablo.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "San Pablo")
remolinoTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Remolino/Remolino.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Remolino")
siekopaiTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Siekopai/Siekopai.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry))
sionaTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Siona/Siona.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Siona")
sinangoeTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Sinangoe/Sinangoe.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Sinangoe")
zabaloTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Zabalo/Zabalo.shp") %>%
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

# plot A: south america
SA <- south_america_sf %>%
    dplyr::filter(continent == "South America") %>%
    ggplot() +
    geom_sf(aes(fill = Ecu), lwd = 0.5) +
    geom_rect(
        aes(
            xmin = e["xmin"], xmax = e["xmax"],
            ymin = e["ymin"], ymax = e["ymax"]
        ),
        fill = NA, color = "#b71212", linewidth = 1.2
    ) +
    theme_void() +
    scale_fill_manual(
        name = "Territory",
        values = c(colors, "A" = "lightyellow", "B" = "white"),
        breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    ) +
    # crop x axis
    coord_sf(
        xlim = c(-80, -35),
        ylim = c(-53, 10)
    )
SA

# plot 2: study area map
studyArea <- south_america_sf %>%
    dplyr::filter(continent == "South America") %>%
    ggplot() +

    # first fill: country background
    geom_sf(aes(fill = Ecu), lwd = 0.5) +
    scale_fill_manual(
        values = c("A" = "lightyellow", "B" = "white"),
        guide = "none"
    ) +
    new_scale_fill() +

    # second fill: natural area raster
    geom_spatraster(data = naturalArea, alpha = 0.2) +
    scale_fill_manual(
        values = c("1" = "forestgreen", "0" = "transparent"),
        labels = c("1" = "Natural Area"),
        breaks = "1",
        name = NULL,
        na.value = "transparent"
    ) +
    new_scale_fill() +

    # third fill: communities
    geom_sf(data = territories_sf, aes(fill = Community), alpha = 0.4, lwd = 0.5) +
    geom_point(data = stations, aes(gps_x, gps_y, fill = Community), alpha = 0.5, pch = 21, size = 1) +
    geom_sf(data = communities_sf, aes(fill = Community), size = 3, pch = 24) +
    geom_text_repel(
        data = communities_sf,
        aes(label = Name, geometry = geometry),
        box.padding = 1.25,
        force = 2,
        force_pull = 0.3,
        min.segment.length = 0,
        max.overlaps = Inf,
        family = "Times",
        stat = "sf_coordinates",
        segment.alpha = 0.4,
        seed = 123
    ) +
    scale_fill_manual(
        name = "Territory",
        values = colors,
        breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    ) +
    scale_y_continuous(position = "right") +
    coord_sf(
        xlim = c(e["xmin"] - 0.05, e["xmax"] + 0.05),
        ylim = c(e["ymin"] - 0.05, e["ymax"] + 0.05),
        label_axes = list(bottom = "E", right = "N")
    ) +
    annotation_scale(
        location = "bl", width_hint = 0.3,
        line_width = 0.5
    ) +
    annotation_north_arrow(
        location = "tr",
        which_north = "true",
        height = unit(1, "cm"),
        width = unit(1, "cm"),
        style = north_arrow_fancy_orienteering
    ) +
    theme_bw() +
    theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
        axis.title = element_blank(),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(family = "Times", colour = "black", size = 16),
        panel.border = element_rect(color = "#b71212", linewidth = 2)
    )
studyArea

# plot C: percent natural area within 25 km
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
natStats <- siteCovariate %>%
    group_by(Community) %>%
    summarize(avgNat = mean(natArea10KM), sdNat = sd(natArea10KM), n = n()) %>%
    mutate(seNat = sdNat / sqrt(n))
natStats$Community <- gsub("Zabalo", "Zábalo", x = natStats$Community)
siteCovariate$Community <- gsub("Zabalo", "Zábalo", x = siteCovariate$Community)

# order the communities by percent of natural cover
natStats <- natStats %>% dplyr::filter(Community %in% orderedCommunities)
natStats$Community <- factor(natStats$Community, levels = orderedCommunities)
natStats

# plot it
natArea <- ggplot(natStats, aes(x = Community, y = avgNat, fill = Community)) +
    geom_bar(stat = "identity") +
    geom_point(
        data = siteCovariate, aes(x = Community, y = natArea10KM),
        position = position_jitter(width = 0.1),
        size = 2, alpha = 0.5, pch = 21
    ) +
    geom_errorbar(aes(ymin = avgNat - sdNat, ymax = avgNat + sdNat), width = 0.2) +

    # put N over each of the bars
    geom_text(aes(label = paste0("N = ", n), y = 0.025), # avgNat + seNat + 0.033
        size = 4, family = "Times", fontface = "bold"
    ) +
    ylab("Percent natural area within 10 km (±SD)") +
    scale_fill_manual(values = colors) +
    ylim(c(0, 1.08)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    theme_bw() +
    theme(
        text = element_text(family = "Times", colour = "black", size = 16),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()
    )
natArea

# violin plot
# 
# p_violin <- ggplot(siteCovariate, aes(x = Community, y = natArea10KM, fill = Community)) +
#     geom_violin(alpha = 0.5, scale = "width", adjust = 1.2) +
#     geom_point(aes(color = Community),
#         position = position_jitter(width = 0.1),
#         size = 0.8, alpha = 1, pch = 21
#     ) +
#     geom_boxplot(
#         width = 0.08, alpha = 0.8,
#         linewidth = 0.3
#     ) +
#     scale_fill_manual(values = colors) +
#     scale_color_manual(values = colors) +
#     labs(x = "Territory", y = "Proportion of natural area within 25km") +
#     theme_bw() +
#     theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Times"),
#         axis.title = element_text(size = 10, family = "Times", face = "bold"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank()
#     ) +
#     ylim(0, 1.05)
# p_violin
# combine plots
# create top panel with SA as inset in top left of studyArea
top_panel <- ggdraw() +
    draw_plot(studyArea, x = 0.09, height = 0.8, scale = 0.83) +
    draw_plot(SA, x = 0.05, y = 0.65, width = 0.3, height = 0.3)
top_panel

# combine top panel with natArea
final_plot <- plot_grid(
    top_panel,
    natArea,
    ncol = 1,
    rel_heights = c(0.7, 0.6),
    labels = c("A", "B"),
    label_size = 18,
    label_fontface = "bold",
    label_fontfamily = "Times"
)
final_plot

# save it
ggsave("Global/Figures/mapPlot.png",
    plot = final_plot,
    width = 8,
    height = 12,
    dpi = 300,
    bg = "white"
)
# pdf
ggsave("Global/Figures/mapPlot.pdf",
    plot = final_plot,
    width = 8,
    height = 12,
    bg = "white"
)
