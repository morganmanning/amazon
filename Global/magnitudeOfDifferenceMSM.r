
############# QUANTIFYING MAGNITUDE OF DIFFERENCE IN MEANS IN OCCUPANCY AT EACH COMMUNITY

setwd("~/Documents/amazon/Global/Data")

# load libraries
require(tidyverse)
require(effectsize)
require(knitr)
require(kableExtra)
require(gridExtra)
require(ggpubr)
require(reshape2)

# load data
load("R Objects/peccaryagoutiMSM_byCommunity.RData")
peccary_agouti <- community_conditional_occupancy
load("R Objects/pacaagoutiMSM_byCommunity.RData")
paca_agouti <- community_conditional_occupancy
load("R Objects/ocelotagoutiMSM_byCommunity.RData")
ocelot_agouti <- community_conditional_occupancy
Traps <- read.csv("AllStationsFormatted.csv")

# one big happy family
bigDF <- rbind(peccary_agouti, paca_agouti, ocelot_agouti)

# want to look per community at the difference in means 
communities <- factor(unique(bigDF$Community),
    levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)
species <- unique(bigDF$Species1)

# N for each community (AKA the number of sites at each community)
nSites <- data.frame(
    Community = communities,
    nSites = c(
        length(unique(Traps$Station[Traps$Community == "Zabalo"])),
        length(unique(Traps$Station[Traps$Community == "Remolino"])),
        length(unique(Traps$Station[Traps$Community == "Sinangoe"])),
        length(unique(Traps$Station[Traps$Community == "San Pablo"])),
        length(unique(Traps$Station[Traps$Community == "Siona"]))
    )
)

# get Hedge's g for each community
HedgesG <- bigDF[,c("Community", "Species1", "Species2")]
HedgesG <- HedgesG[!duplicated(HedgesG), ]
HedgesG$g <- NA
HedgesG$SD1 <- NA
HedgesG$SD2 <- NA
HedgesG$SE1 <- NA
HedgesG$SE2 <- NA

# loop through communities and species
for (i in 1:length(communities)){
    for (j in 1:length(species)){
        # subset data
        subDF <- bigDF[bigDF$Community == communities[i] & bigDF$Species1 == species[j], ]

        # get N
        N <- nSites$nSites[nSites$Community == communities[i]]

        # calculate hedge's g
        SD1 <- subDF$SE[1] * N
        SD2 <- subDF$SE[2] * N
        pooledSD <- sqrt(((N - 1) * (SD1^2) + (N - 1) * (SD2^2)) / (N + N - 2)) # https://www.statisticshowto.com/pooled-standard-deviation/
        #pooledSD <- sqrt((SD1^2 + SD2^2) / 2) # https://www.statisticshowto.com/pooled-standard-deviation/
        mean1 <- subDF$Occupancy[1]
        mean2 <- subDF$Occupancy[2]
        g <- (mean1 - mean2) / pooledSD

        # store results
        HedgesG$g[HedgesG$Community == communities[i] & HedgesG$Species1 == species[j]] <- g
        HedgesG$SD1[HedgesG$Community == communities[i] & HedgesG$Species1 == species[j]] <- SD1
        HedgesG$SD2[HedgesG$Community == communities[i] & HedgesG$Species1 == species[j]] <- SD2
        HedgesG$SE1[HedgesG$Community == communities[i] & HedgesG$Species1 == species[j]] <- subDF$SE[1]
        HedgesG$SE2[HedgesG$Community == communities[i] & HedgesG$Species1 == species[j]] <- subDF$SE[2]
        
        
    }

}

write.csv(HedgesG, "hedges_g.csv")

