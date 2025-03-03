##### PAIRWISE CONDITIONAL PROBABILITY COMPARISON #####
setwd("/Users/morganmanning/Documents/amazon/Global/Data")
setwd("~/Documents/amazon/Global/Data")

################################################################################
# ------------------------------ START UP -------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# load necessary packages
require(dplyr)
require(lubridate)
require(camtrapR)
require(unmarked)
require(ggplot2)
require(rphylopic)
require(knitr)
require(kableExtra)
require(stringr)
require(gridExtra)
require(ggpubr)
require(reshape2)

# load in the necessary data
Data <- read.csv("AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("AllStationsFormatted.csv")
covariates <- read.csv("AllCommunityCovariates.csv")
covariates$Community <- factor(covariates$Community,
    levels = c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
ZABhunting <- read.csv("../../Zabalo/Data/HuntingData2018.csv")

# replace all Mazama species with Mazama sp.
Data$Species <- gsub("Mazama americana", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama nemorivaga", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama gouazoubira", "Mazama sp.", Data$Species)

# get tally of each species at each community
speciesTally <- Data |> 
  group_by(Species, CommunityName) |>
  summarize(nDetections = n()) |>
  #filter(nDetections > 10) |>
  group_by(Species) |>
  mutate(nCommunities = n()) |>
  arrange(desc(nDetections), desc(nCommunities)) |>
  filter(nCommunities == 5) # only pull species that were detected in all four communities

huntingTally <- ZABhunting |> 
  group_by(Species) |>
  summarize(nHunted = n()) |>
  filter(nHunted > 10) |>
  arrange(desc(nHunted))


################################################################################
# ------------------------- DETECTION MATRICES --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# predator/prey; 
  # ocelot (Leopardus pardalis), common opossum (Didelphis marsupialis)
  # include time of day as a covariate (?) for prey since possum is nocturnal, squirrel is diurnal, and ocelot is crepsular

# two most common species (e.g., P(peccary|paca) being lower in disintegrated spaces because people are going to hunt areas with more desirable species; 
  # peccary (Pecari tajacu), paca (Cuniculus paca), black agouti (most spotted but not super hunted, Dasyprocta fuliginosa)

# two species that are going to be competing in the same niche
  # black agouti (Dasyprocta fuliginosa) and green acouchi (Myoprocta pratti)

# for fun
  # (Panthera onca)




# species of interest                       ************* INPUT ***************
species <- c(
            #"Leopardus pardalis",
             #"Didelphis marsupialis", 
             #"Pecari tajacu", 
             "Cuniculus paca", 
             "Dasyprocta fuliginosa" # by FAR the most detected species, but not hunted in ZAB much
             #"Psophia crepitans",
             #"Mazama sp."
             #"Myoprocta pratti", 
             # "Panthera onca" # for funsies!
             ) 
commonNames <- c(
                # "Ocelot", 
                 #"Common opossum", 
                 #"Collared peccary", 
                 "Lowland paca", 
                 "Black agouti"
                 #"Grey-winged trumpeter",
                 #"Brockets"
                 #"Green acouchi",
                 #"Jaguar"
                 ) 
casualNames <- c(
                #"ocelot",
                 #"opossum",
                 #"peccary",
                 "paca",
                 "agouti"
                 #"trumpeter",
                 #"brockets"
                 #"acouchi",
                 #"jaguar"
                 )

# only the interactions we're interested in
interactionsOfInterest <- c(
                            "[ocelot]",
                            "[opossum]",
                            "[peccary]",
                            "[paca]",
                            "[agouti]",
                            "[acouchi]",
                            "[jaguar]",
                            "[ocelot:opossum]",
                            "[agouti:acouchi]",
                            "[peccary:paca]",
                            "[paca:agouti]",
                            "[ocelot:jaguar]"
                            )

# camera operability matrix
Operation <- cameraOperation(CTtable = Traps,
                             stationCol = "Station",
                             cameraCol = "Camera",
                             setupCol = "Setup_date",
                             retrievalCol = "Retrieval_date",
                             hasProblems = TRUE,
                             byCamera = FALSE,
                             allCamsOn = FALSE,
                             camerasIndependent = FALSE,
                             dateFormat = "%Y-%m-%d",
                             writecsv = FALSE)



################################################################################
############################ DETECTION MATRICES ################################
################################################################################
# set up blank lists
detection <- list()

# detection matrices
for (i in 1:length(species)) {
  # occasion length
  occasion = 2 # picked arbitrarily
  # species detection histories for occupancy analyses
  DetHis = detectionHistory(recordTable = Data,
                            camOp = Operation,
                            output = "binary", # binary or count
                            stationCol = "Station",
                            speciesCol = "Species",
                            recordDateTimeCol = "DateTimeOriginal",
                            recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                            day1 = "Station",
                            occasionLength = occasion,
                            datesAsOccasionNames = FALSE,
                            timeZone = "America/Guayaquil",
                            includeEffort = TRUE,
                            scaleEffort = FALSE,
                            #maxNumberDays = 90, #need to think about this
                            species = species[i]) #change species here
  detection[[i]] <- DetHis[["detection_history"]]
  names(detection)[i] <- casualNames[i]
}



################################################################################
########################## MANUAL INPUT REQUIRED ###############################
################################################################################



########################## ALTERNATIVE DIRECTION ###############################
# group species based on heavily hunted vs not really hunted
# this could potentially get at the issue with low detection/overlap

# combine all brocket deer species
# top desirable hunted species based on Zabalo data: peccary (group white-lipped and collared), paca, currasow/guan, agouti
# grouping ideas: all brockets, all birds in family Cracidae (guan/currasow), peccaries, 


wantToGroupSpecies <- "NO" #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (wantToGroupSpecies == "YES") {

    species_groups <- list(
    "Brockets" = c("Mazama americana", "Mazama nemorivaga", "Mazama gouazoubira", "Mazama sp."),
    "Cracidae" = c("Mitu salvini", "Penelope jacquacu"),
    "Peccaries" = c("Pecari tajacu", "Tayassu pecari"),
    "Rodents" = c("Cuniculus paca", "Dasyprocta fuliginosa", "Myoprocta pratti"),
    "Competitors" = c("Leopardus pardalis", "Panthera onca", "Puma concolor")
)
commonNames <- names(species_groups)
casualNames <- names(species_groups)
species <- names(species_groups)

# detection matrices
detection <- list()
for (i in 1:length(commonNames)) {
    group_detection <- list()

    # get a detection matrix for each species within the group
    for (j in 1:length(species_groups[[i]])) {
        # occasion length
        occasion = 2 # picked arbitrarily
        # species detection histories for occupancy analyses
        DetHis = detectionHistory(recordTable = Data,
                                camOp = Operation,
                                output = "binary", # binary or count
                                stationCol = "Station",
                                speciesCol = "Species",
                                recordDateTimeCol = "DateTimeOriginal",
                                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                day1 = "Station",
                                occasionLength = occasion,
                                datesAsOccasionNames = FALSE,
                                timeZone = "America/Guayaquil",
                                includeEffort = TRUE,
                                scaleEffort = FALSE,
                                #maxNumberDays = 90, #need to think about this
                                species = species_groups[[i]][j]) #change species here
        group_detection[[j]] <- DetHis[["detection_history"]]
        names(group_detection)[j] <- species_groups[[i]][j]
    }

    # sum all the group matrices to get one per group with only 1s, 0s, and NAs
    detection[[i]] <- Reduce("+", group_detection)
    names(detection)[i] <- species_groups_names[i]
    detection[[i]] <- ifelse(detection[[i]] > 0, 1, ifelse(detection[[i]] == 0, 0, NA))

}

}
















# clump the detection matrices
clumpedDetections <- list()
for (j in 1:length(detection)) { # for each species
    y <- detection[[j]] # detection history for each species
    clumpEvery <- 2
    nClumpedColumns <- ncol(y) / clumpEvery
    clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = nrow(y))

    clumpStart <- seq(1, ncol(y), by = clumpEvery) # the first column in the clump
    clumpEnd <- seq(clumpEvery, ncol(y), by = clumpEvery) # the last column in the clump

    ### make the clumped matrix
    for (k in 1:nrow(y)) { # for every camera trap station
        for (m in 1:ncol(clumpedMatrix)) {
            if (all(is.na(y[k, clumpStart[m]:clumpEnd[m]])) == TRUE) {
                clumpedMatrix[k, m] <- NA
            } else if (sum(y[k, clumpStart[m]:clumpEnd[m]], na.rm = TRUE) >= 1) {
                clumpedMatrix[k, m] <- 1
            } else {
                clumpedMatrix[k, m] <- 0
            }
        }
    } # the clumped matrix is made
    clumpedDetections[[j]] <- clumpedMatrix
    names(clumpedDetections)[j] <- casualNames[j]
}



################################################################################
########################## MANUAL INPUT REQUIRED ###############################
################################################################################

# DO YOU WANT TO PROCEED WITH THE CLUMPED DETECTION MATRICES?
clumped <- "NO" # "YES" or "NO"
if (clumped == "YES") {
    detection <- clumpedDetections
}

################################################################################
################################################################################
################################################################################



# identify empty rows (i.e. stations with no detections)
empty_rows <- vector()
for (i in 1:length(detection)){
  allNArow <- c(as.numeric(which(rowSums(is.na(detection[[i]]))==ncol(detection[[i]]))))
  empty_rows <- append(empty_rows, allNArow)
}

# remove empty rows if they're empty in all species 
detectionWithoutBlanks <- list()
if(all(table(empty_rows) == length(species)) & 
all(covariates$Station == rownames(detection[[1]]))) { # if every species has the same empty rows
  rowsToCull <- unique(empty_rows)
  for (i in 1:length(detection)){
    detectionWithoutBlanks[[i]] <- detection[[i]][-rowsToCull,]
  }
  covariatesWithoutBlanks <- covariates[-rowsToCull,]
  names(detectionWithoutBlanks) <- names(detection)
}
covariatesWithoutBlanks$Year <- as.factor(covariatesWithoutBlanks$Year)
covariatesWithoutBlanks$Station <- as.factor(covariatesWithoutBlanks$Station)
covariatesWithoutBlanks$Community <- as.factor(covariatesWithoutBlanks$Community)
covariatesWithoutBlanks$NearestCommunity <- as.factor(covariatesWithoutBlanks$NearestCommunity)

# make the unmarked frame
umf <- unmarkedFrameOccuMulti(
    y = detectionWithoutBlanks,
    siteCovs = covariatesWithoutBlanks, 
    maxOrder = 2 # max number of species interactions
)
summary(umf)
str(umf)









################################################################################
############################ OCCUPANCY MODELING ################################
################################################################################

# ALL COMMUNITIES TOGETHER:
# Percentage natural area (see formattingLULC.R for buffer size)
# Community as a covariate
# Average monthly rainfall from July-November (kg/m^2/s)
# Temperature (C)
# Distance to a water source (m)

match_detection <- c("Community")
match_occupancy <- c("RainfallScaled", "percentNatural", "DistToWater", "TemperatureScaled", "Year")
# excluded "Community" from occupancy covariates bc it correlated with percentNatural (per chisq.test())


############# BEST OCCUPANCY MODEL PER SPECIES USING THE BEST DETECTION MODEL
# every possible combination of variables
combos <- sapply( seq(length(match_occupancy)), function(i) {
    as.list(as.data.frame(combn( x = match_occupancy, m = i)))
})
combos <- unlist(combos, recursive=FALSE)
  
# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))
forms <- as.vector(c(forms, "~ 1")) # add the null model
occupancyFormulas <- forms
  
############### COMBINE DETECTION WITH ALL POSSIBLE OCCUPANCY PREDICTORS AND MODEL
occupancyModelsList <- list()
for (j in 1:length(species)) {
df <- data.frame(detection = '~ Community',
                occupancy = occupancyFormulas)
occupancyModelsList[[j]] <- c(paste(df$detection, df$occupancy, sep = " "))
}

# run occupancy unmarked model for all models
allModels <- list()
ufoList <- list()
for (j in 1:length(species)) { # for all the species
    occupancyMods <- list()
    ufo <- unmarkedFrameOccu(detectionWithoutBlanks[[j]],
        siteCovs = covariatesWithoutBlanks,
        obsCovs = NULL
    )
    ufoList[[j]] <- ufo # store all the unmarked frames
    for (m in 1:length(occupancyModelsList[[j]])) {
        test <- occu(formula(occupancyModelsList[[j]][m]), ufo)
        occupancyMods[[m]] <- occu(formula(occupancyModelsList[[j]][m]), ufo,
            control = 10000,
            starts = c(rep(0, length(test@opt$par)))
        )
    }
    names(occupancyMods) <- 1:length(occupancyMods)
    allModels[[j]] <- occupancyMods
    print(paste0("Finishing species ", j, " out of ", length(species)))
}

# Make a data frame to show the model names and their AICs
modelAICs <- list() # all models and their AICs
topModels <- list() # only models within 2 AIC of lowest AIC
for (j in 1:length(allModels)) { # for each species
    df <- data.frame(
        ModelName = NA,
        AIC = NA,
        diffFromBest = NA
    )
    for (m in 1:length(allModels[[j]])) {
        df[m, 1] <- as.character(c(allModels[[j]][[m]]@formula))
        df[m, 2] <- allModels[[j]][[m]]@AIC
    }
    df <- df[order(df$AIC), ]
    df$diffFromBest <- df$AIC - min(df$AIC)
    modelAICs[[j]] <- df

    # only the best
    ANTM <- subset(df, diffFromBest <= 2)
    topModels[[j]] <- ANTM
}
names(topModels) <- casualNames

# create a new dataframe with list titles as a new column
topModels_df <- data.frame(
    Species = character(),
    ModelName = character(),
    AIC = numeric(),
    diffFromBest = numeric(),
    stringsAsFactors = FALSE
)

for (j in 1:length(topModels)) {
    speciesX <- casualNames[j]
    models <- topModels[[j]]
    n <- nrow(models)

    # create a temporary dataframe for each species
    temp_df <- data.frame(
        Species = rep(speciesX, n),
        ModelName = models$ModelName,
        AIC = models$AIC,
        diffFromBest = models$diffFromBest
    )

    # append the temporary dataframe to the main dataframe
    topModels_df <- rbind(topModels_df, temp_df)
}

# save the table
kbl(topModels_df) %>%
  kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
  save_kable(file = "../Figures/MultispeciesModeling/speciesBestModels.png", zoom = 2)


# plot raw detection
for (j in 1:length(ufoList)) { # for every species
    ufo <- ufoList[[j]]
    colnames(ufo@y) <- 1:ncol(ufo@y)
    meltedComb <- melt(ufo@y)
    meltedComb$value <- as.factor(meltedComb$value)

    # plot clumped
    p <- ggplot(meltedComb, aes(Var2, Var1, fill = value)) +
        geom_tile(colour = "gray50") +
        scale_fill_manual(values = c("gray75", "red3"), na.value = "white", name = "") +
        scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
        scale_alpha_identity(guide = "none") +
        coord_equal(expand = 0) +
        xlab("Time (1 unit = ~2 days)") +
        ylab("Camera trap site") +
        ggtitle(commonNames[j]) +
        theme_bw() +
        theme(plot.title = element_text(size = 25, hjust = 0.5))
    ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[j], "RawDetection.png"),
        plot = p, width = 7, height = 8
    )
}

# for running very simple model with just brockets and agouti
if (length(species) == 2) {
    # both species 
    ufo1 <- ufoList[[1]]
    ufo2 <- ufoList[[2]]

    # add them together
    combinedDetection <- ufo1@y + ufo2@y
    colnames(combinedDetection) <- 1:ncol(combinedDetection)
    meltedComb <- melt(combinedDetection)
    meltedComb$value <- as.factor(meltedComb$value)

    # plot clumped
    p <- ggplot(meltedComb, aes(Var2, Var1, fill = value)) +
        geom_tile(colour = "gray50") +
        scale_fill_manual(values = c("gray75","#f77c7c", "red3"), na.value = "white", name = "") +
        scale_x_continuous(breaks = seq(0, ncol(ufo2@y), by = 2)) +
        scale_alpha_identity(guide = "none") +
        coord_equal(expand = 0) +
        xlab("Time (1 unit = ~2 days)") +
        ylab("Camera trap site") +
        ggtitle(paste0("Overlapping detection of ", casualNames[1], " and ", casualNames[2])) +
        theme_bw() +
        theme(plot.title = element_text(size = 25, hjust = 0.5))
    ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[1], casualNames[2], "OverlappingRawDetection.png"),
        plot = p, width = 7, height = 8
    )

}

################################################################################
################################ MODELING ######################################
################################################################################
# look at occupancy covariates
topModels_df$OccupancyFormula <- gsub("~Community", "", topModels_df$ModelName)
topModels_df[c("Species", "OccupancyFormula")]

# extract the words/letter strings
modelCovariates <- strsplit(gsub("[ ~]", "", topModels_df$OccupancyFormula), "[ +~]")
names(modelCovariates) <- topModels_df$Species

# make modelCovariates into a dataframe with the list titles as a new column
modelCovariates_df <- data.frame(Species = character(),
                 Covariates = character())

for (i in 1:length(modelCovariates)) {  
  # create a temporary dataframe for each species
  temp_df <- data.frame(Species = rep(names(modelCovariates)[i], length(modelCovariates[[i]])),
            Covariates = modelCovariates[[i]])
  
  # append the temporary dataframe to the main dataframe
  modelCovariates_df <- rbind(modelCovariates_df, temp_df)
}

# pull the best covariates for each species and make a formula
bestCovariates <- modelCovariates_df |> 
  unique() |>
  subset(Covariates  != "1") |>
  mutate(ID = row_number()) |>
  group_by(Species) |>
  summarize(BestCovariates = paste("~", paste(Covariates, collapse = "+"))) 
bestCovariates <- bestCovariates[match(casualNames, bestCovariates$Species), ]
bestCovariates

# verify that the species are in identical orders
if (all(bestCovariates$Species == names(detectionWithoutBlanks))) {
    
    # use the formulas to run occuMulti (with ~Community as the detection formula)
    best_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            bestCovariates$BestCovariates,
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    # following recommendations of https://groups.google.com/g/unmarked/c/0gSJXk_Ew94/m/Xto_7YnTBAAJ
    best_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            bestCovariates$BestCovariates,
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        control = list(maxit = 20000),
        #method = "Nelder-Mead",
        starts = rep(0, length(best_multispecies_model@opt$par)),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    print("Finished optimal model :)")

    # the null model
    null_multispecies_model <- occuMulti(
        detformulas = rep("~ 1", length(species)),
        stateformulas = c(
            rep("~1", (length(bestCovariates$BestCovariates))),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    null_multispecies_model <- occuMulti(
        detformulas = rep("~ 1", length(species)),
        stateformulas = c(
            rep("~1", (length(bestCovariates$BestCovariates))),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2,
        control = list(maxit = 20000),
        #method = "Nelder-Mead",
        starts = rep(0, length(null_multispecies_model@opt$par))
    )
    print("Finished null model :)")

    # just percent natural as a covariate since it was in all the best models
    natural_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ percentNatural", length(bestCovariates$BestCovariates)),
            rep("~ percentNatural", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    natural_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ percentNatural", length(bestCovariates$BestCovariates)),
            rep("~ percentNatural", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2,
        control = list(maxit = 20000),
        #method = "Nelder-Mead",
        starts = rep(0, length(natural_multispecies_model@opt$par))
    )
    print("Finished percent natural area-only model :)")

    # just community as a covariate since it was in all the best models
    community_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ Community", length(bestCovariates$BestCovariates)),
            rep("~ Community", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    community_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ Community", length(bestCovariates$BestCovariates)),
            rep("~ Community", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2,
        control = list(maxit = 20000),
        # method = "Nelder-Mead",
        starts = rep(0, length(community_multispecies_model@opt$par))
    )
    print("Finished community-only model :)")


    # global model
    global_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ RainfallScaled + percentNatural + DistToWater + TemperatureScaled + Year", length(bestCovariates$BestCovariates)),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
    global_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ RainfallScaled + percentNatural + DistToWater + TemperatureScaled + Year", length(bestCovariates$BestCovariates)),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2,
        control = list(maxit = 20000),
        #method = "Nelder-Mead",
        starts = rep(0, length(global_multispecies_model@opt$par))
    )
    print("Finished global model :)")

}
# at this point, the best covariates were used for each species' occupancy with "Community" as the detection covariate 
# and no covariates for interactions

summary(null_multispecies_model)
summary(best_multispecies_model)
summary(global_multispecies_model)
summary(natural_multispecies_model)
summary(community_multispecies_model)
# it's not looking good, girl...


# add penalty to model: Penalized likelihood estimation
set.seed(123)
best_mod_penalty <- optimizePenalty(best_multispecies_model, penalties = c(0.5, 1))
summary(best_mod_penalty)
natural_mod_penalty <- optimizePenalty(natural_multispecies_model, penalties = c(0.5, 1))
summary(natural_mod_penalty)
community_mod_penalty <- optimizePenalty(community_multispecies_model, penalties = c(0.5, 1))
summary(community_mod_penalty)


save(null_multispecies_model, best_multispecies_model, global_multispecies_model,
    natural_multispecies_model, community_multispecies_model, 
    best_mod_penalty, natural_mod_penalty, community_mod_penalty, 
    casualNames, commonNames, species, umf,
    file = "R Objects/multispeciesModels.RData"
)



