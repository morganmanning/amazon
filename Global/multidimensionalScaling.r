# Goal: do the multidimensional scaling

setwd("~/amazon/Global/Data")

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
library(knitr)
library(kableExtra)

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
keepColumns <- c(
    "TerritoryArea", "distanceToNearestCommunity",
    "airTemp", "airTempSD", "rainfall", "rainfallSD", "MeanElevation", "MeanDistToWater"
)
excludeColumns <- c(
    "X", "Y", "OperatingDays", "MeanTemperature", "MeanDistToWater",
    "DaysHuntingPerMonthDry", "DaysHuntingPerMonthWet", "DaysFishingPerMonthWet",
    "DaysFishingPerMonthDry", "PercentPopWhoHunt", "PercentPopWhoFish", 
    "humiditySD", "airTempSD", "rainfallSD", "rootMoistureSD",
    # new removals
    "humidity", "rootMoisture", "shannonIndex", "simpsonIndex", "nSpecies", "nIndiv"
)
communityCovariatesRemoved <- communityCovariates[, (names(communityCovariates) %in% keepColumns)]
communityCovariatesRemoved <- scale(communityCovariatesRemoved)

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
NMDS <- metaMDS(distance_matrix, trymax = 50)

# add arrows showing the direction of the covariates
ef <- envfit(NMDS, as.data.frame(communityCovariatesRemoved), na.rm = TRUE)
ef
labelFactor <- 2.5
arrows1 <- ef$vectors$arrows |> as_tibble(rownames = "community")
as_tibble(NMDS$points, rownames = "community") |>
    ggplot(aes(x = MDS1, y = MDS2, label = community)) +
    #geom_point()
    xlim(min(NMDS$points[, 1]) * 2, max(NMDS$points[, 1]) * 1.2) +
    geom_tile(
        aes(fill = community, color = "black", alpha = 0.2),
        width = 1,
        height = 0.2,
        col = "black"
    ) +
    scale_fill_manual(values = colors) +
    geom_text(aes(fontface = "bold"), size = 7) +
    geom_segment(data = arrows1, aes(x = labelFactor * NMDS1, y = labelFactor * NMDS2, xend = 0, yend = 0), arrow = arrow(ends = "first")) +
    geom_label_repel(data = arrows1, aes(x = labelFactor * NMDS1, y = labelFactor * NMDS2), nudge_y = .1, segment.colour = NA) +
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

# make a table of the envfit results
ef_df <- as.data.frame(ef$vectors$arrows) %>%
    mutate(
        r2 = ef$vectors$r,
        pval = ef$vectors$pvals,
        Variable = rownames(ef$vectors$arrows)
    ) %>%
    select(NMDS1, NMDS2, r2, pval) %>%
    mutate(
        sig = case_when(
            pval <= 0.001 ~ "***",
            pval <= 0.01 ~ "**",
            pval <= 0.05 ~ "*",
            #pval <= 0.1 ~ ".",
            TRUE ~ ""
        ),
        pval_fmt = sprintf("%.3f%s", pval, sig),
        r2_fmt = sprintf("%.3f", r2)
    ) %>%
    select(NMDS1, NMDS2, r2_fmt, pval_fmt)

# make the table
kbl(
    ef_df,
    col.names = c("NMDS1", "NMDS2", expression(r^2), "p-value"),
    escape = FALSE,
    format = "html"
) %>%
    kable_classic(full_width = TRUE, html_font = "Times New Roman") %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, bold = TRUE) %>%   # bolds the header row
    save_kable(
        "../Figures/MultispeciesModeling/NMDS_envfit_table.png",
        zoom = 1.5
    )
