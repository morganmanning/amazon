# Goal: do the multidimensional scaling

setwd("~/Documents/amazon/Global/Data")

# source of elevation data: https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d
# source for rainfall, temperature, humidity, root moisture: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary

################################################################################
################################################################################
################################################################################

# load in necessary packages
require(dplyr)

# load in necessary data
communityCovariates <- read.csv("CommunityLevelCovariates.csv")
rownames(communityCovariates) <- communityCovariates$Community
communityCovariates$Community <- NULL



################################################################################
################################################################################
######################### MULTIDIMENSIONAL SCALING #############################
################################################################################
################################################################################





