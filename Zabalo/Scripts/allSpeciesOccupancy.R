############################################################################
######################## OCCUPANCY OF ALL SPECIES ##########################
############################################################################
# requirements: input files must be detection matrices (e.g., "Deer.csv" in Zabalo)

setwd("~/Documents/amazon/Zabalo/Data")

require(unmarked)
require(dplyr)
require(ggplot2)
require(reshape2)


############################################################################
############################# LOAD DATA ####################################
############################################################################


  ##############################################
  ##### !!! insert all .csv file names !!! #####
  ##############################################
speciesDFs <- c(
  "CollaredPeccary.csv",
  "Deer.csv",
  "Paca.csv"
)

# read in all .csv files and remove any ID rows or columns 
dfList <- list()
for (i in 1:length(speciesDFs)) {
  species <- read.csv(speciesDFs[i])
  if (any(species[,1] > 1, na.rm = TRUE) | 
      any(species[1,] > 1, na.rm = TRUE)) { # if first row or column is ID
    if (any(species[1,] > 1, na.rm = TRUE)) { # row
      print(paste("ERROR: Suspected ID *row* removed from", speciesDFs[i]))
    } else if (any(species[,1] > 1, na.rm = TRUE)) { # column
      species <- data.frame(species[,-1], row.names = species[,1]) # row names = stations
    } else next 
  } else { # if first row or column is NOT ID
    print(paste("ERROR: No ID rows or columns removed from", speciesDFs[i]))
  }
  dfList[[i]] <- species
}

# name dfs after their associated species
listTitles <- gsub('.csv', '', speciesDFs)
names(dfList) <- listTitles

# Zabalo occupancy covariates
load('R Objects/siteCovs2018.RData')
names(siteCovariate)

# stations and hunting info
stations <- read.csv("Stations2018.csv")
hunting <- read.csv("HuntingData2018.csv")
cameraRecords <- read.csv("RecordTable2018.csv")


############################################################################
######################## FORMAT FOR UNMARKED ###############################
############################################################################

# make all the data frames in the list into matrices
for (i in 1:length(dfList)) {
  dfList[[i]] <- as.matrix(dfList[[i]])
  dfList[[i]] <- (dfList[[i]][ order(as.numeric(row.names(dfList[[i]]))), ]) #order matters
}

# clump all matrices according to their best clumping factor
source("../Scripts/optimalClumping.R")

clumpedMatrixList <- list()
for (i in 1:length(dfList)) {
  y <- dfList[[i]]
  clumpEvery <- as.numeric(best_clumping_factor(y))
  nClumpedColumns <- ncol(y)/clumpEvery
  clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = nrow(y)) 
  
  clumpStart <- seq(1, ncol(y), by = clumpEvery) # the first column in the clump
  clumpEnd <- seq(clumpEvery, ncol(y), by = clumpEvery) # the last column in the clump
  
  ### make the clumped matrix
  for (k in 1:30){
    for (j in 1:ncol(clumpedMatrix)){
      if(all(is.na(y[k, clumpStart[j]:clumpEnd[j]])) == TRUE) {
        clumpedMatrix[k,j] <- NA
      } else if (sum(y[k, clumpStart[j]:clumpEnd[j]], na.rm=TRUE) >= 1) {
        clumpedMatrix[k,j] <- 1
      } else {
        clumpedMatrix[k,j] <- 0
      }
    }
  } # the clumped matrix is made
  clumpedMatrixList[[i]] <- clumpedMatrix
}
names(clumpedMatrixList) <- listTitles


# make unmarked df for each species
ufoList <- list()
for (i in 1:length(clumpedMatrixList)){
  clumpedMatrix <- clumpedMatrixList[[i]]
  ufo <- unmarkedFrameOccu(clumpedMatrix, 
                           siteCovs = siteCovariate,
                           obsCovs = NULL)
  ufoList[[i]] <- ufo
}
names(ufoList) <- listTitles

# make unmarked df for each species without clumping
unclumpedUFOList <- list()
for (i in 1:length(dfList)){
  unclumpedMatrix <- dfList[[i]]
  ufo <- unmarkedFrameOccu(unclumpedMatrix, 
                           siteCovs = siteCovariate,
                           obsCovs = NULL)
  unclumpedUFOList[[i]] <- ufo
}
names(unclumpedUFOList) <- listTitles

# make detection plots
# clumped
for (i in 1:length(ufoList)) {
  ufo <- ufoList[[i]]
  colnames(ufo@y) <- 1:ncol(ufo@y)
  meltedComb <- melt(ufo@y)
  meltedComb$value <- as.factor(meltedComb$value)
  
  #plot
  ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
    geom_tile(colour = "gray50") +
    scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
    scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    xlab(paste("Time (1 unit = ~", 
               as.numeric(best_clumping_factor(dfList[[i]]))*2,
               # number of columns/time steps divided by the clumping factor, times 2
               "days)")) +
    ylab("Camera trap site") +
    theme_bw()
  ggsave(filename = gsub(" ", "", paste("clumped_", names(ufoList)[i], ".png")), 
         width = 8, height = 4, 
         path = "../Figures/DetectionPlots")

}

# unclumped
for (i in 1:length(unclumpedUFOList)) {
  ufo <- unclumpedUFOList[[i]]
  colnames(ufo@y) <- 1:ncol(ufo@y)
  meltedComb <- melt(ufo@y)
  meltedComb$value <- as.factor(meltedComb$value)
  
  #plot
  ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
    geom_tile(colour = "gray50") +
    scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
    scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    xlab("Time (1 unit = 2 days)") +
    ylab("Camera trap site") +
    theme_bw()
  ggsave(filename = gsub(" ", "", paste("rawDetection_", names(unclumpedUFOList)[i], ".png")), 
         width = 8, height = 4, 
         path = "../Figures/DetectionPlots")
  
}


############################################################################
############################## DETECTION ###################################
############################################################################

# vector of variables to consider for detection
match_variables <- c("Effort", "Habitat", "OnTrail")

# every possible combination of variables
combos <- sapply( seq(length(match_variables)), function(i) {
  as.list(as.data.frame(combn( x = match_variables, m = i)))
})
combos <- unlist(combos, recursive=FALSE)

# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))

detectionFormulas <- as.vector(c(forms, "~ 1"))

# get the best detection formula and stick with that for the occupancy formulas
tempDF <- data.frame(detection = detectionFormulas,
                     occupancy = "~ 1")
allDetectionFormulas <- paste(tempDF$detection, tempDF$occupancy, sep = " ")

# see which detection model is the best
bestDetectionModels <- data.frame(species = listTitles,
                                  bestDetectionModel = NA)
for (j in 1:length(ufoList)) {
  detectionMods <- list()
  temp <- data.frame(Model = allDetectionFormulas,
                     AIC = NA)
  for(i in 1:length(allDetectionFormulas)) {
    test <- occu(formula(allDetectionFormulas[[i]]), ufoList[[j]])
    detectionMods[[i]] <- occu(formula(allDetectionFormulas[[i]]), ufoList[[j]], 
                               control = 10000, 
                               starts = c(rep(0, length(test@opt$par)))) 
    temp$AIC[i] <- detectionMods[[i]]@AIC
  }
  # select the model with the lowest AIC
  bestDetectionModels$bestDetectionModel[j] <- detectionFormulas[which(temp$AIC == min(temp$AIC))]
  }



############################################################################
############################## OCCUPANCY ###################################
############################################################################

# vector of variables 
match_variables <- c("Habitat", "Trail.Distance", "HuntingIntensity")

# every possible combination of variables
combos <- sapply( seq(length(match_variables)), function(i) {
  as.list(as.data.frame(combn( x = match_variables, m = i)))
})
combos <- unlist(combos, recursive=FALSE)

# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))

forms <- as.vector(c(forms, "~ 1"))

occupancyFormulas <- forms


############################################################################
############################### MODELS #####################################
############################################################################
occupancyModelsList <- list()
for (i in 1:length(ufoList)) {
  df <- data.frame(detection = bestDetectionModels$bestDetectionModel[i],
                                         occupancy = occupancyFormulas)
  occupancyModelsList[[i]] <- c(paste(df$detection, df$occupancy, sep = " "))
  }

# run occupancy unmarked model for all models
allModels <- list()
for (j in 1:length(ufoList)) { # for all the species
  occupancyMods <- list()
  for(i in 1:length(occupancyModelsList[[j]])) {
    test <- occu(formula(occupancyModelsList[[j]][i]), ufoList[[j]])
    occupancyMods[[i]] <- occu(formula(occupancyModelsList[[j]][i]), ufoList[[j]], 
                               control = 10000, 
                               starts = c(rep(0, length(test@opt$par)))) 
  }
  names(occupancyMods) <- 1:length(occupancyMods)
  allModels[[j]] <- occupancyMods
}

save(allModels, file = 'R Objects/AllSpeciesModels.RData')

# remove all models with missing SE/z/p-value or that didn't converge
noMissingMods <- list()
for(j in 1:length(allModels)) {
  eachSpecies <- list()
  for(i in 1:length(allModels[[j]])){
    
    modSum <- summary(allModels[[j]][[i]])
    
    if(anyNA(modSum$state$SE)==FALSE & allModels[[j]][[i]]@opt$convergence != 1){
      eachSpecies[[i]] <- allModels[[j]][[i]]
    } else {
      next
    }
  }
  eachSpecies[sapply(eachSpecies, is.null)] <- NULL
  noMissingMods[[j]] <- eachSpecies
}

# remove all models with ridiculous standard errors
noWackyMods <- list()
for(j in 1:length(noMissingMods)) {
  eachSpecies <- list()
  for(i in 1:length(noMissingMods[[j]])){
    
    modSum <- summary(noMissingMods[[j]][[i]])
    
    if((abs(modSum$state$SE[1]) < modSum$state$Estimate[1]) == TRUE){
      eachSpecies[[i]] <- noMissingMods[[j]][[i]]
    } else {
      next
    }
  }
  eachSpecies[sapply(eachSpecies, is.null)] <- NULL
  noWackyMods[[j]] <- eachSpecies
}
names(noWackyMods) <- listTitles
save(noWackyMods, file = 'R Objects/noWackyMods.RData')


############################################################################
########################### MODEL SELECTION ################################
############################################################################

# Make a data frame to show the model names and their AICs
modelAICs <- list()
topModels <- list() # only models within 2 AIC of lowest AIC
for(j in 1:length(noWackyMods)) {
  df <- data.frame(ModelName = NA,
                   AIC = NA,
                   diffFromBest = NA)
  for(i in 1:length(noWackyMods[[j]])){
    df[i,1]<- as.character(c(noWackyMods[[j]][[i]]@formula))
    df[i,2]<- noWackyMods[[j]][[i]]@AIC
  }
  df <- df[order(df$AIC),]
  df$diffFromBest <- df$AIC - min(df$AIC)
  modelAICs[[j]] <- df
  
  # only the best
  ANTM <- subset(df, diffFromBest <= 2)
  topModels[[j]] <- ANTM
}
names(topModels) <- listTitles

# take these best models and put them into a list
bestModsFitLists <- list()
for (j in 1:length(topModels)){
  speciesBestMods <- list()
  for(i in 1:nrow(topModels[[j]])) {
    test <- occu(formula(topModels[[j]]$ModelName[i]), ufoList[[j]])
    
    speciesBestMods[[i]] <- occu(formula(topModels[[j]]$ModelName[i]), ufoList[[j]], 
                          control = 10000, 
                          starts = c(rep(0, length(test@opt$par)))) 
    #print(speciesBestMods[[i]])
  } 
  bestModsFitLists[[j]] <- fitList(speciesBestMods)
}
names(bestModsFitLists) <- listTitles


############################################################################
##################### PLOTTING MODEL PREDICTIONS ###########################
############################################################################


############### TRAIL DISTANCE ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
trailDF <- data.frame(
                 Habitat = as.factor("Upland"), # most common habitat
                 Trail.Distance = seq(min(siteCovariate$Trail.Distance), 
                                      max(siteCovariate$Trail.Distance), 
                                      length.out = 100),
                 HuntingIntensity = mean(siteCovariate$HuntingIntensity)
) 

# predict
trailPredictions <- list()
for (j in 1:length(bestModsFitLists)) {
  mods <- bestModsFitLists[[j]]
  unmarkedPred <- predict(mods, type = 'state', new = trailDF, appendData = TRUE)
  predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                    lower = unmarkedPred$lower,
                                    upper = unmarkedPred$upper,
                                    trailDF)
  trailPredictions[[j]] <- predictionDataFrame
}
names(trailPredictions) <- listTitles
save(trailPredictions, file = 'R Objects/trailPredictions.RData')


############### HUNTING INTENSITY ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
huntingDF <- data.frame(
                 Habitat = as.factor("Upland"), # most common habitat
                 Trail.Distance = mean(siteCovariate$Trail.Distance),
                 HuntingIntensity = seq(min(siteCovariate$HuntingIntensity), 
                                        max(siteCovariate$HuntingIntensity), 
                                        length.out = 100)
) 

# predict
huntingPredictions <- list()
for (j in 1:length(bestModsFitLists)) {
  mods <- bestModsFitLists[[j]]
  unmarkedPred <- predict(mods, type = 'state', new = huntingDF, appendData = TRUE)
  predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                    lower = unmarkedPred$lower,
                                    upper = unmarkedPred$upper,
                                    huntingDF)
  huntingPredictions[[j]] <- predictionDataFrame
}
names(huntingPredictions) <- listTitles
save(huntingPredictions, file = 'R Objects/huntingPredictions.RData')



###### average occupancy and detection across sites
# Model-averaged prediction of occupancy and confidence interval

df <- data.frame(
  Habitat = siteCovariate$Habitat,
  HuntingIntensity = siteCovariate$HuntingIntensity,
  Trail.Distance = siteCovariate$Trail.Distance
) 

estimatedParameters <- list()
for (j in 1:length(bestModsFitLists)) {
  # state = occupancy
  unmarkedPredOcc <- predict(bestModsFitLists[[j]], type = 'state', new = df, appendData = TRUE)
  # det = detection
  unmarkedPredDet <- predict(bestModsFitLists[[j]], type = 'det', new = df, appendData = TRUE)
  estimatedParameters[[j]] <- data.frame(avgOccupancy = mean(unmarkedPredOcc$Predicted),
                                         avgOccupancySE = mean(unmarkedPredOcc$SE),
                                         occupancyRange = range(unmarkedPredOcc$Predicted),
                                         avgDetection = mean(unmarkedPredDet$Predicted),
                                         avgDetectionSE = mean(unmarkedPredDet$SE),
                                         detectionRange = range(unmarkedPredDet$Predicted))
}
names(estimatedParameters) <- listTitles



# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)
mean(unmarkedPred$Predicted); mean(unmarkedPred$SE)
range(unmarkedPred$Predicted)

unmarkedPred <- predict(bestMods, type = 'det', new = df, appendData = TRUE)
mean(unmarkedPred$Predicted); mean(unmarkedPred$SE)
range(unmarkedPred$Predicted)
