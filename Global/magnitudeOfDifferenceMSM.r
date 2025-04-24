
############# QUANTIFYING MAGNITUDE OF DIFFERENCE IN MEANS IN OCCUPANCY AT EACH COMMUNITY


# load data
load("R Objects/peccary_agoutiMSM_byCommunity.RData")
peccary_agouti <- community_conditional_occupancy
load("R Objects/paca_agoutiMSM_byCommunity.RData")
paca_agouti <- community_conditional_occupancy
load("R Objects/ocelot_agoutiMSM_byCommunity.RData")
ocelot_agouti <- community_conditional_occupancy


# want to look per community at the difference in means 
