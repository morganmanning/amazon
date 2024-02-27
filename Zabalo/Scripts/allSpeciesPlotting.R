
############################################################################
#################### PLOTTING OCCUPANCY OF ALL SPECIES #####################
############################################################################
# requirements: must have run allSpeciesOccupancy.R to get necessary R objects

setwd("~/Documents/amazon/Zabalo/Data")
library(rphylopic)
require(ggpubr)
require(dplyr)
require(ggplot2)


############################################################################
############################# LOAD DATA ####################################
############################################################################

load('R Objects/trailPredictions.RData')
load('R Objects/huntingPredictions.RData')
load('R Objects/speciesNames.RData')

# get limits of covariates for each species (for plotting)
xLimits <- data.frame(species = speciesNames,
                      trailMinX = NA,
                      trailMaxX = NA,
                      huntingMinX = NA,
                      huntingMaxX = NA)
for (j in 1:length(speciesNames)) {
  speciesTrailPrediction <- trailPredictions[[j]]
  xLimits[j,"trailMinX"] <- min(speciesTrailPrediction$Trail.Distance)
  xLimits[j,"trailMaxX"] <- max(speciesTrailPrediction$Trail.Distance)
  
  speciesHuntingPrediction <- huntingPredictions[[j]]
  xLimits[j,"huntingMinX"] <- min(speciesHuntingPrediction$HuntingIntensity)
  xLimits[j,"huntingMaxX"] <- max(speciesHuntingPrediction$HuntingIntensity)

}


# make a plot for hunting intensity and trail distance for each species j
for (j in 1:length(speciesNames)) {
  
  # Distance to a trail prediction
  speciesTrailPrediction <- trailPredictions[[j]]
  ggplot(speciesTrailPrediction, aes(x = Trail.Distance, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'green4', alpha = 0.5, linetype = "dashed") +
    geom_path(linewidth = 1) +
    labs(x = "Distance to a trail (scaled)", y = "Occupancy probability") +
    ggtitle(label = speciesNames[j]) +
    theme_classic() +
    coord_cartesian(ylim = c(0,1), xlim = c(xLimits[j,"trailMinX"], xLimits[j,"trailMaxX"])) +
    theme(text = element_text(family = "HelveticaNeue", colour = "black"),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, face = 'italic'))
  
  ggsave(filename = gsub(" ", "", paste("trailOccupancy_", speciesNames[j], ".png")), 
         width = 8, height = 4, 
         path = "../Figures/OccupancyPrediction")
  
  # Hunting intensity
  speciesHuntingPrediction <- huntingPredictions[[j]]
  ggplot(speciesHuntingPrediction, aes(x = HuntingIntensity, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
    geom_path(linewidth = 1) +
    labs(x = "Hunting intensity (scaled)", y = "Occupancy probability") +
    ggtitle(label = speciesNames[j]) +
    theme_classic() +
    coord_cartesian(ylim = c(0,1), xlim = c(xLimits[j,"huntingMinX"], xLimits[j,"huntingMaxX"])) +
    theme(text = element_text(family = "HelveticaNeue", colour = "black"),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, face = 'italic'))
  ggsave(filename = gsub(" ", "", paste("huntingOccupancy_", speciesNames[j], ".png")), 
         width = 8, height = 4, 
         path = "../Figures/OccupancyPrediction")
}




plot(1:10, 1:10)


############# NOTES
# Get a single image uuid for a species
uuid <- get_uuid(name = "Cuniculus paca", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(pacaHunting <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.25))

# figure panel
arranged <- ggarrange(grid::nullGrob(), pacaRiver, pacaTrail,
                      #labels = c("  ", "B", "C"),
                      ncol = 3, nrow = 1)
annotate_figure(arranged, top = text_grob("Lowland paca", face = "bold", size = 16))








