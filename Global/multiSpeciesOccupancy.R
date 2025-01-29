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

# load in the necessary data
Data <- read.csv("AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("AllStationsFormatted.csv")
covariates <- read.csv("AllCommunityCovariates.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
ZABhunting <- read.csv("../../Zabalo/Data/HuntingData2018.csv")

# get tally of each species at each community
speciesTally <- Data |> 
  group_by(Species, CommunityName) |>
  summarize(nDetections = n()) |>
  #filter(nDetections > 10) |>
  group_by(Species) |>
  mutate(nCommunities = n()) # |>
  # filter(nCommunities == 4) # only pull species that were detected in all four communities

huntingTally <- ZABhunting |> 
  group_by(Species) |>
  summarize(nHunted = n()) |>
  filter(nHunted > 10)


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
species <- c("Leopardus pardalis",
             "Didelphis marsupialis", 
             "Pecari tajacu", 
             "Cuniculus paca", 
             "Dasyprocta fuliginosa", # by FAR the most detected species, but not hunted in ZAB much
             "Myoprocta pratti", 
             "Panthera onca") # for funsies!
commonNames <- c("Ocelot", 
                 "Common opossum", 
                 "Collared peccary", 
                 "Lowland paca", 
                 "Black agouti",
                 "Green acouchi",
                 "Jaguar") # listTitles
casualNames <- c("ocelot",
                 "opossum",
                 "peccary",
                 "paca",
                 "agouti",
                 "acouchi",
                 "jaguar")

# only the interactions we're interested in
interactionsOfInterest <- c("[ocelot]",
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
                            "[ocelot:jaguar]")

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
  occasion = 10 # picked arbitrarily
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
str(detectionWithoutBlanks)
nrow(covariatesWithoutBlanks)

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
match_occupancy <- c("Rainfall", "percentNatural", "DistToWater", "Temperature")
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
for (j in 1:length(species)) { # for all the species
    occupancyMods <- list()
    ufo <- unmarkedFrameOccu(detectionWithoutBlanks[[j]],
        siteCovs = covariatesWithoutBlanks,
        obsCovs = NULL
    )
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
  save_kable(file = "../Figures/MultispeciesModeling/7speciesBestModels.png", zoom = 2)

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
    multispecies_model <- occuMulti(
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

    # just percent natural as a covariate since it was in all the best models
    natural_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ percentNatural", length(bestCovariates$BestCovariates)),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )

    # global model
    global_multispecies_model <- occuMulti(
        detformulas = rep("~ Community", length(species)),
        stateformulas = c(
            rep("~ Rainfall + percentNatural + DistToWater + Temperature", length(bestCovariates$BestCovariates)),
            rep("~ 1", ((length(species)^2 + length(species)) / 2) - length(species))
        ),
        # (n^2+n)/2 is the addition version of a factorial
        # use null for all animal interactions
        data = umf,
        maxOrder = 2
    )
}
# at this point, the best covariates were used for each species' occupancy with "Community" as the detection covariate 
# and no covariates for interactions

summary(null_multispecies_model)
# it's not looking good, girl...







################################################################################
########################## OCCUPANCY PREDICTIONS ###############################
################################################################################
# make a dataframe with predictions for all species
all_predictions <- data.frame()

for (i in 1:length(casualNames)){
    # get predictions for each species
    preds <- predict(null_multispecies_model, type = "state", species = "paca")
    all_predictions <- rbind(predictions, preds[1,])
}
all_predictions$Species <- casualNames

# plot null model predictions for occupancy
plot(1:length(species), all_predictions$Predicted,
    #ylim = c(0.1, 0.4),
    #xlim = c(0.5, 3.5), 
    pch = 19, cex = 1.5, xaxt = "n",
    xlab = "", ylab = "Marginal occupancy and 95% CI"
)
axis(1, at = 1:length(species), labels = all_predictions$Species)

# CIs
top <- 0.1
for (i in 1:length(species)) {
    segments(i, all_predictions$lower[i], i, all_predictions$upper[i])
    segments(i - top, all_predictions$lower[i], i + top)
    segments(i - top, all_predictions$upper[i], i + top)
}
