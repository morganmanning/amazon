
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

# one big happy family
bigDF <- rbind(peccary_agouti, paca_agouti, ocelot_agouti)

# want to look per community at the difference in means 
communities <- factor(unique(bigDF$Community),
    levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)
species <- unique(bigDF$Species1)

for (i in 1:length(communities)){
    for (j in 1:length(species)){
        # subset data
        subDF <- bigDF[bigDF$Community == communities[i] & bigDF$Species1 == species[j], ]
        
    }

}

write.csv(peccary_agouti, "peccary_agouti.csv")

