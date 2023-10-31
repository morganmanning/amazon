
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


############################################################################
############################### PECCARY ####################################
############################################################################

### not edited
# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = Trail.Distance, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'green4', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to a trail (scaled)", y = "Occupancy probability") +
  #ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 2.05)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryTrail <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))



#### not edited
# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = HuntingIntensity, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Hunting intensity (scaled)", y = "Occupancy probability") +
  # ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 3)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryHunting <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))

# figure panel
arranged <- ggarrange(peccaryComm, peccaryRiver, peccaryTrail, peccaryHunting,
                      #labels = c("A", "B", "C"),
                      ncol = 4, nrow = 1)
annotate_figure(arranged, top = text_grob("Collared peccary", face = "bold", size = 16))




############################################################################
############################### BROCKET ####################################
############################################################################

############### HUNTING INTENSITY ##################

# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = HuntingIntensity, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Hunting intensity (scaled)", y = "Occupancy probability") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 3)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Mazama pandora", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(deerHunting <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.15, ysize = 0.25))


# figure panel
arranged <- ggarrange(deerComm, deerRiver,
                      #labels = c("A", "B"),
                      ncol = 3, nrow = 1)
annotate_figure(arranged, top = text_grob("Brown brocket", face = "bold", size = 16))






############################################################################
######################@@######### PACA #####################################
############################################################################

############### TRAIL DISTANCE ##################

# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = Trail.Distance, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'green4', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to a trail (scaled)", y = "Occupancy probability") +
  #ggtitle(label = "Lowland paca") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 2.05)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Cuniculus paca", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(pacaTrail <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.75))



############### HUNTING INTENSITY ##################

# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = HuntingIntensity, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Hunting intensity (scaled)", y = "Occupancy probability") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 3)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

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








