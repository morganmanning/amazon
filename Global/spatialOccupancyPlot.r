
# GOAL: plot occupancy by community to visualize spatial patterns

################################################################################
############################ LOAD AND PREPARE DATA #############################
################################################################################

# prep
rm(list = ls())
setwd("~/Documents/amazon")

# load libraries
library(tidyverse)
library(sf)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)
library(terra)
library(MASS)
library(unmarked)
require(tidyterra) # for geom_spatraster
require(ggnewscale)

# load master data
load("Global/Data/R Objects/masterBestModsFitLists.RData")
# CONTAINS:
# commonNames ("Ocelot", "Black agouti", etc)
# speciesNames (scientific names)
# masterBestModsFitLists[[1]][[j]] where j is the species
# masterBestModsOutputs[[1]][[j]][[m]] where j is the species and m is the model number
# siteCovariate (dataframe with site covariates, including community names and site names)
# stations (dataframe with station info)
# modelAveragesDF_clean (dataframe with model averaged estimates for all species and covariates)

# load stations/covariate data
stations <- read.csv("Global/Data/AllStationsFormatted.csv")
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
agriculture <- rast("Global/Data/Community-level Covariates/DEM/agriculture.tif")
agriculture <- as.factor(agriculture)


# format stations data
stations <- stations %>%
    dplyr::select(c(Station, gps_x, gps_y, CommunityName)) %>%
    distinct()
stations$CommunityName <- gsub("Zabalo", "Zábalo", x = stations$CommunityName)

# order communities by natural cover
orderedCommunities <- c(siteCovariate %>%
    group_by(Community) %>%
    summarize(perc = mean(PercentNaturalScaled)) %>%
    arrange(desc(perc)) %>%
    dplyr::select(Community))$Community
orderedCommunities <- gsub("Zabalo", "Zábalo", x = orderedCommunities)

stations <- stations %>% dplyr::filter(CommunityName %in% orderedCommunities)
stations$Community <- factor(stations$CommunityName, levels = orderedCommunities)
stations$CommunityName <- NULL

# define colors
colors <- c(
    "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
)

# load South America shapefile
south_america_sf <- ne_countries(
    scale = "large",
    continent = "south america",
    returnclass = "sf"
)

# load territory shapefiles
sanPabloTerritory <- st_read("Global/Data/Community-level Covariates/Territories/SanPablo/SanPablo.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "San Pablo")
remolinoTerritory <- st_read("Global/Data/Community-level Covariates/Territories/Remolino/Remolino.shp") %>%
    dplyr::select(c(Name, Shape_Leng, Shape_Area, geometry)) %>%
    dplyr::mutate(Community = "Remolino")
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
    sionaTerritory, sinangoeTerritory, zabaloTerritory
)

# load community points
remSP <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Secoya.kml")
remSP <- remSP %>%
    dplyr::mutate(Community = ifelse(Name == "San Pablo" | Name == "Waiya" | Name == "Bellavista", "San Pablo", "Remolino"))
siona <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Siona.kml") %>%
    dplyr::mutate(Community = "Siona")
sinangoe <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Sinangoe.kml") %>%
    dplyr::mutate(Community = "Sinangoe")
zabalo <- sf::st_read("Global/Data/Community-level Covariates/Community Points/Zabalo.kml") %>%
    dplyr::mutate(Community = "Zábalo") %>%
    dplyr::mutate(Name = ifelse(Name == "Cofan Zabalo", "Zábalo", Name))

communities_sf <- rbind(remSP, siona, sinangoe, zabalo)

# bounding box
coords <- as.matrix(stations[, c("gps_x", "gps_y")])
e <- as.vector(ext(coords))
e["xmin"] <- e["xmin"] - 0.1
e["ymin"] <- e["ymin"] - 0.5
e["xmax"] <- e["xmax"] + 0.1
e["ymax"] <- e["ymax"] + 0.5


################################################################################
########################## OCCUPANCY PREDICTION ################################
################################################################################

# species names
species_list <- speciesNames

# prepare newdata for prediction
newdata <- siteCovariate %>%
    dplyr::select(
        Station, Community, RainfallScaled, NatArea20KMScaled,
        TemperatureScaled, DistToComm
    ) %>%
    distinct()

all_occupancy_predictions <- list()
# loop through each species
for (i in 1:length(species_list)) {
    sp_name <- species_list[i]

    # list of best model fits for this species
    best_models <- masterBestModsFitLists[[1]][[i]]

    if (is.null(best_models) || length(best_models) == 0) {
        message(paste("No models for", sp_name))
        next
    }

    # store predictions from each model
    model_predictions <- list()

    for (j in 1:length(best_models)) {
        tryCatch(
            {
                # get predictions from this model
                pred <- predict(best_models, type = "state", newdata = newdata)

                # extract predicted occupancy
                model_predictions[[j]] <- pred$Predicted
            },
            error = function(e) {
                message(paste("Error predicting for", sp_name, "model", j, ":", e$message))
            }
        )
    }

    # average predictions across models
    if (length(model_predictions) > 0) {
        pred_matrix <- do.call(cbind, model_predictions)
        avg_occupancy <- rowMeans(pred_matrix, na.rm = TRUE)

        # results
        result_df <- newdata %>%
            dplyr::select(Station, Community) %>%
            mutate(
                species = sp_name,
                occupancy = avg_occupancy
            )

        all_occupancy_predictions[[sp_name]] <- result_df
    }
}

# combine all predictions
occupancy_predictions <- bind_rows(all_occupancy_predictions)
occupancy_predictions$Community <- gsub("Zabalo", "Zábalo", x = occupancy_predictions$Community)

# merge with station coordinates
occupancy_spatial <- occupancy_predictions %>%
    left_join(stations, by = c("Station", "Community")) %>%
    dplyr::filter(!is.na(gps_x) & !is.na(gps_y)) # remove NA coordinates

################################################################################
######################### PLOT OCCUPANCY PER SITE ##############################
################################################################################

for (i in 1:length(unique(occupancy_spatial$species))) {
    sp <- unique(occupancy_spatial$species)[i]

    species_data <- occupancy_spatial %>%
        dplyr::filter(species == sp)

    # create base map
    ecuador_sf <- south_america_sf %>%
        dplyr::filter(name_en == "Ecuador")

    # extract covariates from best models for this species
    best_models <- masterBestModsOutputs[[1]][[i]]

    all_covariates <- c()

    if (!is.null(best_models) && length(best_models) > 0) {
        for (j in 1:length(best_models)) {
            tryCatch(
                {
                    # get formula from the model
                    model_formula <- best_models[[j]]@formula
                    state_formula <- model_formula[[3]] # occupancy component
                    formula_terms <- all.vars(state_formula)
                    formula_terms <- formula_terms[formula_terms != "~"]
                    all_covariates <- c(all_covariates, formula_terms)
                },
                error = function(e) {
                    message(paste("Could not extract formula from model", j, "for", sp))
                }
            )
        }
    }

    # get unique covariates
    unique_covariates <- unique(all_covariates)

    # readable covariate names
    covariate_map <- c(
        "RainfallScaled" = "Rainfall",
        "NatArea20KMScaled" = "Natural Area (within 20km)",
        "TemperatureScaled" = "Temperature",
        "DistToComm" = "Distance to Community",
        "Community" = "Community"
    )

    # map to readable names
    readable_covariates <- sapply(unique_covariates, function(x) {
        if (x %in% names(covariate_map)) {
            return(covariate_map[x])
        } else {
            return(x)
        }
    })

    # list out covariates
    if (length(readable_covariates) > 0) {
        caption_text <- paste("Covariates:", paste(readable_covariates, collapse = ", "))
    } else {
        caption_text <- "Null model used"
    }

    # change title based on whether community was a covariate
    if (sd(species_data$occupancy, na.rm = TRUE) == 0) {
        plotTitle <- paste("Predicted Model-averaged Occupancy:", sp)
        plotSubtitle <- "(Community was not a covariate in the best models)"
    } else {
        plotTitle <- paste("Predicted Model-averaged Occupancy:", sp)
        plotSubtitle <- NULL
    }

    # plot it
    p_points <- ggplot() +
        geom_sf(data = ecuador_sf, fill = "lightyellow", color = "gray50", lwd = 0.5) +
        geom_sf(
            data = territories_sf, aes(color = Community),
            fill = NA, alpha = 0.7, lwd = 1
        ) +
        geom_spatraster(data = agriculture, alpha = 0.25) +
        scale_fill_manual(
            values = c("1" = "#773d2b", "0" = "transparent"),
            labels = c("1" = "Agriculture"),
            breaks = "1",
            name = NULL,
            na.value = "transparent"
        ) +
        new_scale_fill() +
        geom_point(
            data = species_data,
            aes(x = gps_x, y = gps_y, fill = occupancy),
            pch = 21, alpha = 0.8, color = "black"
        ) +
        scale_fill_gradientn(
            colors = c("darkblue", "lightblue", "yellow", "orange", "red"),
            values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
            limits = c(0, 1),
            name = "Occupancy\nProbability"
        ) +
        new_scale_fill() +
        geom_sf(data = communities_sf, aes(fill = Community), size = 3, pch = 24) +
        scale_fill_manual(
            values = colors,
            name = "Community",
            guide = "none"
        ) +
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
        scale_color_manual(
            values = colors,
            name = "Territory"
        ) +
        scale_size_continuous(range = c(2, 10), guide = "none") +
        coord_sf(
            xlim = c(e["xmin"] - 0.05, e["xmax"] + 0.05),
            ylim = c(e["ymin"] - 0.05, e["ymax"] + 0.05)
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
        labs(
            title = plotTitle,
            subtitle = plotSubtitle,
            caption = caption_text
        ) +
        theme_bw() +
        theme(
            panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
            panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
            axis.title = element_blank(),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5),
            plot.caption = element_text(size = 8, hjust = 0, face = "italic"),
            text = element_text(family = "Times", colour = "black"),
            legend.position = "right"
        )

    # save plot
    ggsave(paste0("Global/Figures/SingleSpeciesModeling/occupancy_points_", gsub(" ", "_", sp), ".png"),
        p_points,
        width = 10, height = 8, dpi = 300
    )

    print(paste("Saved point plot for", sp))
}



