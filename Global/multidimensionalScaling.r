# Goal: do the multidimensional scaling

setwd("~/Documents/amazon/Global/Data")

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d
# source for rainfall, temperature, humidity, root moisture: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary

################################################################################
################################################################################
################################################################################

# load in necessary packages
require(dplyr)
require(vegan)
require(ggplot2)
require(ggrepel)

# load in necessary data
communityCovariates <- read.csv("CommunityLevelCovariates.csv")
rownames(communityCovariates) <- gsub("Zabalo", "Zábalo", communityCovariates$Community)
communityCovariates$Community <- NULL



################################################################################
################################################################################
######################### MULTIDIMENSIONAL SCALING #############################
################################################################################
################################################################################

# Columns to exclude
excludeColumns <- c(
    "X", "Y", "OperatingDays", "MeanTemperature", "MeanDistToWater",
    "DaysHuntingPerMonthDry", "DaysHuntingPerMonthWet", "DaysFishingPerMonthWet",
    "DaysFishingPerMonthDry", "PercentPopWhoHunt", "PercentPopWhoFish", 
    "humiditySD", "airTempSD", "rainfallSD", "rootMoistureSD"
)
communityCovariatesRemoved <- communityCovariates[, !(names(communityCovariates) %in% excludeColumns)]


# Calculate the distance matrix
distance_matrix <- dist(communityCovariatesRemoved)

# Perform MDS analysis
mds <- cmdscale(distance_matrix)
plot(mds, type = "n")
text(mds, labels = rownames(communityCovariates))

# data frame for ggplot
x <- as.data.frame(mds)$V1
y <- -as.data.frame(mds)$V2 # reflect so North is at the top
mds_df <- data.frame(x = x, y = y, label = rownames(mds))

# plot using ggplot
colors <- c(
    "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
)
ggplot(mds_df, aes(x = x, y = y, label = label)) +
    geom_point() +
    geom_text(aes(label = label), vjust = -0.5) + # add labels
    theme_bw() +
    theme(aspect.ratio = 1) 




# this is how I found a way to add the arrows to the plot
####### the 'vegan' way based on https://andrewirwin.github.io/data-viz-notes/lessons/122-mds.html
NMDS <- metaMDS(distance_matrix, trace = 0)

# add arrows showing the direction of the covariates
ef <- envfit(NMDS, communityCovariatesRemoved, na.rm = TRUE)
ef

arrows1 <- ef$vectors$arrows |> as_tibble(rownames = "community")
as_tibble(NMDS$points, rownames = "community") |>
    ggplot(aes(x = MDS1, y = MDS2, label = community)) +
    xlim(min(NMDS$points[, 1]) * 1.3, max(NMDS$points[, 1]) * 1.3) +
    geom_tile(
        aes(fill = community, color = "black", alpha = 0.5),
        width = 210000000,
        col = "black"
    ) +
    scale_fill_manual(values = colors) +
    geom_text(aes(fontface = "bold"), size = 5) +
    geom_segment(data = arrows1, aes(x = 50000000 * NMDS1, y = 50000000 * NMDS2, xend = 0, yend = 0)) +
    geom_label_repel(data = arrows1, aes(x = 50000000 * NMDS1, y = 50000000 * NMDS2),
     max.overlaps = 15,
     min.segment.length = .5, box.padding = 0.9) +
    theme_bw() +
    theme(aspect.ratio = 1,
        legend.position = "none")
ggsave(
    filename = paste0("../Figures/MultispeciesModeling/multiDimensionalScaling.png"),
    width = 12, height = 12
)

# base plot for comparison
plot(NMDS, type = "text", cex = 1.5)
plot(ef)
