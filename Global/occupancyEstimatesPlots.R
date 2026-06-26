
################################################################################
############# OCCUPANCY ESTIMATES ACROSS COMMUNITIES AND SPECIES ###############
################################################################################
# mission: plot estimates for each species in each community
rm(list = ls())
setwd("~/Documents/amazon")

# load packages
require(unmarked)
require(dplyr)
require(ggplot2)
require(reshape2)
# require(rphylopic)
require(ggpubr)
require(sf)
require(spData)
#require(ggmagnify)
require(terra)
require(mapdata)
require(kableExtra)
require(knitr)
require(lubridate)
#require(tictoc)
#require(MuMIn)
require(tidyverse)
require(ggimage)
require(magick)

#tic() # time it

# input
communities <- "Global"
communitiesAbrv <- "All"
#speciesNames <- c("Pecari tajacu", "Mazama americana", "Cuniculus paca", "Psophia crepitans")
#commonNames <- c("Collared peccary", "Red brocket", "Lowland paca", "Grey-winged trumpeter") # listTitles
#speciesNames <- c("Eira barbara")
speciesNames <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Leopardus pardalis")
#commonNames <- c("Tayra") 
commonNames <- c("Collared peccary", "Brockets", "Lowland paca", "Grey-winged trumpeter", "Black agouti", "Nine-banded armadillo", "Great tinamou", "Ocelot") 
  # paca = Cuniculus paca
  # brocket = Mazama americana
  # collared peccary = Pecari tajacu 
  # trumpeter = Psophia crepitans
  # brown four-eyed possum = Metachirus nudicaudatus (#1 species in SGE)
  # black agouti = Dasyprocta fuliginosa (#2 species in SGE)





################################################################################
################################################################################
######################### QUESTIONS BEFORE PROCEEDING ##########################
################################################################################
################################################################################


# DO YOU WANT TO PLOT BASED ON THE GROUPED OR UNGROUPED DATA?
proceedWithClumpedData <- "YES"

# DO YOU WANT TO FORCE THE GROUPING TO BE EVERY FOUR DAYS FOR ALL SPECIES?
force4DayGrouping <- "NO"

# DO YOU WANT TO USE THE GLOBAL MODELS FOR EACH SPECIES?
useGlobalModels <- "NO"


################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
## Check correlation between covariates ##

# load site covariates for all communities
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
siteCovariate$Rainfall <- siteCovariate$Rainfall*1000 # convert to grams/m^2/s  

# plot correlation matrix
correlation <- siteCovariate %>%
  select(RainfallScaled, ag10KM, natArea10KM, DistToWater, TemperatureScaled, DistToComm)
correlationMatrix <- cor(correlation, use = "complete.obs")
correlationMatrix
# nat area and ag are highly negatively correlated (-0.96), so proceeded with just nat area




################################################################################
########################## CONDENSING TIME STEPS ###############################
################################################################################
# mission: make unmarked occupancy frames (ufo) and 
# plot the detection across clumped vs. regular time steps


# output lists
clumpPlotMasterList <- list()
unclumpPlotMasterList <- list()
ufoMasterList <- list()
unclumpedUFOMasterList <- list()

# find optimal time step clumping and plot it
for (i in 1:length(communities)) {
  # make a list to put all detection histories into
  detHistory <- list()
  
    for (j in 1:length(speciesNames)){
      csv <- paste0(communitiesAbrv[i], gsub(" ", "", speciesNames[j]), ".csv")
      species <- read.csv(paste0(communities[i], "/Data/", csv))
      species <- species[,-1] # remove the column of ID names
      
      # save each detection history in a list
      detHistory[[j]] <- species    
      }
 
  # now, with the list of species' detection histories, proceed for each community
  names(detHistory) <- commonNames
  
  # stations info
  stations <- read.csv(paste0(communities[i], "/Data/", communitiesAbrv[i], "StationsFormatted.csv"))
  cameraRecords <- read.csv(paste0(communities[i], "/Data/", communitiesAbrv[i], "IndependentRecordsFormatted.csv"))
  # replace all Mazama species with Mazama sp.
  cameraRecords$Species <- gsub("Mazama americana", "Mazama sp.", cameraRecords$Species)
  cameraRecords$Species <- gsub("Mazama nemorivaga", "Mazama sp.", cameraRecords$Species)
  cameraRecords$Species <- gsub("Mazama gouazoubira", "Mazama sp.", cameraRecords$Species) # replace all Mazama species with Mazama sp.

  
  # import site covariates for each community
  if (communities[i] == "Sinangoe"){
    siteCovariate <- data.frame(DistToComm = c(scale(stations$Community/1000))) # scale
      } else if (communities[i] == "Zabalo") {
        load('Zabalo/Data/R Objects/siteCovs2018.RData') # loads 'siteCovariate'
        } else if (communities[i] == "Global") {
         siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
         siteCovariate$Rainfall <- siteCovariate$Rainfall*1000 # convert to grams/m^2/s
         siteCovariate$Community <- factor(siteCovariate$Community, 
                                            levels = c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona"))
         
        } else {
          siteCovariate <- NULL # no covariates for remaining communities 
        }
  
  # turn all species detection histories for community i into matrices
  for (j in 1:length(detHistory)) {
    detHistory[[j]] <- detHistory[[j]][ order(as.numeric(row.names(detHistory[[j]]))), ] #order matters
    detHistory[[j]] <- as.matrix(detHistory[[j]])
  }
  
  # clump all matrices according to their best clumping factor
  source("./Zabalo/Scripts/optimalClumping.R")
  
  clumpedMatrixList <- list()
  nDaysGroupedPerSpecies <- data.frame(
    Species = speciesNames,
    CommonNames = commonNames,
    DaysGrouped = NA
  )
  for (j in 1:length(detHistory)) { # for each species
    y <- detHistory[[j]] # detection history for each species

    # if you want to force grouping by 2: 
    if (force4DayGrouping == "YES") {
      clumpEvery <- 4
    } else {
      clumpEvery <- as.numeric(best_clumping_factor(detHistory[[j]], maximumClumpingFactor = 20)[1])
    }

    # store what each species was grouped by
    nDaysGroupedPerSpecies$DaysGrouped[j] <- clumpEvery

    # make a list of clumped matrices
    nClumpedColumns <- ncol(y)/clumpEvery
    clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = nrow(y)) 
    
    clumpStart <- seq(1, ncol(y), by = clumpEvery) # the first column in the clump
    clumpEnd <- seq(clumpEvery, ncol(y), by = clumpEvery) # the last column in the clump
    
    ### make the clumped matrix
    for (k in 1:nrow(y)){ # for every camera trap station
      for (m in 1:ncol(clumpedMatrix)){
        if(all(is.na(y[k, clumpStart[m]:clumpEnd[m]])) == TRUE) {
          clumpedMatrix[k,m] <- NA
        } else if (sum(y[k, clumpStart[m]:clumpEnd[m]], na.rm=TRUE) >= 1) {
          clumpedMatrix[k,m] <- 1
        } else {
          clumpedMatrix[k,m] <- 0
        }
      }
    } # the clumped matrix is made
    clumpedMatrixList[[j]] <- clumpedMatrix
  }
  names(clumpedMatrixList) <- commonNames # now we have a matrix of clumped detection histories
  
  # make unmarked df for each species with clumping
  ufoList <- list()
  for (j in 1:length(clumpedMatrixList)){
    clumpedMatrix <- clumpedMatrixList[[j]]
    ufo <- unmarkedFrameOccu(clumpedMatrix, 
                             siteCovs = siteCovariate,
                             obsCovs = NULL)
    ufoList[[j]] <- ufo
    }
  names(ufoList) <- commonNames
  
  # make unmarked df for each species without clumping
  unclumpedUFOList <- list()
  for (j in 1:length(detHistory)){
    unclumpedMatrix <- detHistory[[j]]
    ufo <- unmarkedFrameOccu(unclumpedMatrix, 
                             siteCovs = siteCovariate,
                             obsCovs = NULL)
    unclumpedUFOList[[j]] <- ufo
  }
  names(unclumpedUFOList) <- commonNames
  
  ## make detection plots
  # clumped
  clumpPlotList <- list()
  for (j in 1:length(ufoList)) { # for every species
    ufo <- ufoList[[j]]
    colnames(ufo@y) <- 1:ncol(ufo@y)
    meltedComb <- melt(ufo@y)
    meltedComb$value <- as.factor(meltedComb$value)
    
    # plot clumped
    clumpPlotList[[j]] <- ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
      geom_tile(colour = "gray50") +
      scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
      scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
      scale_alpha_identity(guide = "none") +
      coord_equal(expand = 0) +
      xlab(paste("Time (1 unit = ~", 
                 as.numeric(best_clumping_factor(detHistory[[j]]))*2,
                 # number of columns/time steps divided by the clumping factor, times 2
                 "days)")) +
      ylab("Camera trap site") +
      ggtitle(paste(communities[i], commonNames[j])) +
      theme_bw() +
      theme(plot.title = element_text(size = 25, hjust = 0.5))
  }
  
  # unclumped
  unclumpPlotList <- list()
  for (j in 1:length(unclumpedUFOList)) {
    ufo <- unclumpedUFOList[[j]]
    colnames(ufo@y) <- 1:ncol(ufo@y)
    meltedComb <- melt(ufo@y)
    meltedComb$value <- as.factor(meltedComb$value)
    
    #plot
    unclumpPlotList[[j]] <- ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
      geom_tile(colour = "gray50") +
      scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
      scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
      scale_alpha_identity(guide = "none") +
      coord_equal(expand = 0) +
      xlab("Time (1 unit = 2 days)") +
      ylab("Camera trap site") +
      ggtitle(paste(communities[i], commonNames[j])) +
      theme_bw() +
      theme(plot.title = element_text(size = 25, hjust = 0.5))
   
  }
  
  # output lists
  clumpPlotMasterList[[i]] <- clumpPlotList
  unclumpPlotMasterList[[i]] <- unclumpPlotList
  ufoMasterList[[i]] <- ufoList
  unclumpedUFOMasterList[[i]] <- unclumpedUFOList 
  
}
names(clumpPlotMasterList) <- communities
names(unclumpPlotMasterList) <- communities
names(ufoMasterList) <- communities
names(unclumpedUFOMasterList) <- communities









######################## CLUMPING VERSUS RAW DETECTION #########################
################################################################################
################################################################################
################################################################################

# BASED ON ANSWER AT THE BEGINNING OF THE SCRIPT
print(paste0("You said ", proceedWithClumpedData, " to proceeding with clumped data"))
print(paste0("You said ", force4DayGrouping, " to forcing 2-day grouping for all species"))

if (proceedWithClumpedData == "NO") {
    ufoMasterList <- unclumpedUFOMasterList
}




################################################################################
################################################################################
################################################################################










###############################################################################
############################ OCCUPANCY MODELING ################################
################################################################################

# ZABALO:
  # detection: c("Effort", "Habitat")
  # occupancy: c("Habitat", "Trail.Distance", "HuntingIntensity")

# SINANGOE:
  # detection: c("1")
  # occupancy: c("DistToComm")

# SIONA/SIEKOPAI:
  # detection: c("1")
  # occupancy: c("1)

# ALL COMMUNITIES TOGETHER:
  # Percentage natural area (see formattingLULC.R for buffer size)
  # Community as a covariate
  # Average monthly rainfall from July-November (kg/m^2/s)
  # Temperature (C)
  # Distance to a water source (m)
  # Distance to the community (m)
  # Distance to a community (km)

# output lists:
masterBestofTheBest <- list() # occu output for the best model for each species
masterTopModels <- list() # df of the models within 2 AIC of the top model
masterBestModsFitLists <- list() # fit list of the best models for model averaging
masterBestModsOutputs <- list()
masterBestModsOutputs <- list()

for (i in 1:length(communities)) {
  
  # set the variables to be considered for America's Next Top Model
  if (communities[i] == "Zabalo") {
    match_detection <- c("Effort", "Habitat")
    match_occupancy <- c("Habitat", "Trail.Distance", "HuntingIntensity")
  } else if (communities[i] == "Sinangoe"){
    match_detection <- c("1")
    match_occupancy <- c("DistToComm")
  } else if (communities[i] == "Siona" | communities[i] == "Siekopai"){
    match_detection <- c("1")
    match_occupancy <- c("1")
  } else if (communities[i] == "Global") {
    match_detection <- c("Community", "DaysEffortScaled")
    match_occupancy <- c("Community", "RainfallScaled", "NatArea10KMScaled",
                         "DistToWater", "TemperatureScaled", "DistToComm")
  }
  
  
  ############### BEST DETECTION MODEL PER SPECIES
  # every possible combination of variables
  combos <- sapply(seq(length(match_detection)), function(k) {
    as.list(as.data.frame(combn(x = match_detection, m = k)))
    })
  combos <- unlist(combos, recursive=FALSE)
  
  # all combinations of variables into formulas
  forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))
  detectionFormulas <- as.vector(c(forms, "~ 1"))
  
  # get the best detection formula and stick with that for the occupancy formulas
  tempDF <- data.frame(detection = detectionFormulas,
                       occupancy = "~ 1")
  allDetectionFormulas <- paste(tempDF$detection, tempDF$occupancy, sep = " ")
  bestDetectionModels <- data.frame(species = commonNames,
                                    bestDetectionModel = NA)
  
  # see which detection model is the best with null occupancy
  for (j in 1:length(speciesNames)){
    # make a df off all detection formulas and their AIC for each species
    detectionMods <- list()
    temp <- data.frame(Model = allDetectionFormulas,
                       AIC = NA)
    
    # run occu model for every detection formula per species per community to compare AIC
    for(m in 1:length(allDetectionFormulas)) {
      test <- occu(formula(allDetectionFormulas[[m]]), ufoMasterList[[i]][[j]])
      detectionMods[[m]] <- occu(formula(allDetectionFormulas[[m]]), ufoMasterList[[i]][[j]], 
                                 control = 10000, 
                                 starts = c(rep(0, length(test@opt$par)))) 
      temp$AIC[m] <- detectionMods[[m]]@AIC
    }
    # select the model with the lowest AIC
    bestDetectionModels$bestDetectionModel[j] <- detectionFormulas[which(temp$AIC == min(temp$AIC))]
    
    # over ride and use null if species is ocelot
    if (commonNames[j] == "Ocelot") {
      bestDetectionModels$bestDetectionModel[j] <- "~ 1"
    }
    # output is dataframe for community i with the best detection formulas for each species
      
  }
  
  
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
  
  
  
  ############### COMBINE BEST DETECTION WITH ALL POSSIBLE OCCUPANCY PREDICTORS AND MODEL
  occupancyModelsList <- list()
  for (j in 1:length(speciesNames)) {
    df <- data.frame(detection = bestDetectionModels$bestDetectionModel[j],
                     occupancy = occupancyFormulas)
    occupancyModelsList[[j]] <- c(paste(df$detection, df$occupancy, sep = " "))
  }
  
  # run occupancy unmarked model for all models
  allModels <- list()
  modelAICs <- list() # all models and their AICs
  topModels <- list() # only models within 2 AIC of lowest AIC
  for (j in 1:length(speciesNames)) { # for all the species
    # if ocelot, force null model
    if (commonNames[j] == "Ocelot") {
      occupancyModelsList[[j]] <- paste("~ 1", "~ 1")
    }
    occupancyMods <- list()
    for(m in 1:length(occupancyModelsList[[j]])) {
      test <- occu(formula(occupancyModelsList[[j]][m]), ufoMasterList[[i]][[j]])
      occupancyMods[[m]] <- occu(formula(occupancyModelsList[[j]][m]), ufoMasterList[[i]][[j]], 
                                 control = 10000, 
                                 starts = c(rep(0, length(test@opt$par)))) 
    }
    names(occupancyMods) <- 1:length(occupancyMods)
    allModels[[j]] <- occupancyMods

    # # make a df taking the model name and AIC
    # # but automatically put an NA if running the model outputs an error
    # df <- data.frame(ModelName = rep(NA, length(occupancyModelsList[[j]])),
    #                  AIC = NA,
    #                  diffFromBest = NA)
    
    # # if (is.na(summary(allModels[[j]][[m]])$state$SE[1]) == TRUE) {
    # #     df[m, 1] <- NA
    # #     df[m, 2] <- NA
    # # } else {
    # #     df[m, 1] <- as.character(c(allModels[[j]][[m]]@formula))
    # #     df[m, 2] <- allModels[[j]][[m]]@AIC
    # # }

    # # if the model outputs an error, put NA (use tryCatch)
    # tryCatch(
    #     suppressWarnings({
    #         model <- allModels[[j]][[m]]
    #         model_summary <- summary(model)

    #         if (is.na(summary(allModels[[j]][[m]])$state$SE[1]) == TRUE) {
    #             df[m, 1] <- NA
    #             df[m, 2] <- NA
    #         } else {
    #             df[m, 1] <- as.character(c(allModels[[j]][[m]]@formula))
    #             df[m, 2] <- allModels[[j]][[m]]@AIC
    #         }
    #     }),
    #     error = function(e) {
    #         df[m, 1:2] <- NA
    #     }
    # )


    # # order the DF by AIC
    # df <- df[order(df$AIC), ]
    # df <- df[!is.na(df$ModelName), ]
    # df$diffFromBest <- df$AIC - min(df$AIC)
    # modelAICs[[j]] <- df

    # # only the best
    # ANTM <- subset(df, diffFromBest <= 2)
    # topModels[[j]] <- ANTM

    # print a message at the end of each species
    print(paste("Done with", j, "out of", length(speciesNames), "species :)"))
  }
  
  # Make a data frame to show the model names and their AICs
  modelAICs <- list() # all models and their AICs
  topModels <- list() # only models within 2 AIC of lowest AIC
  for(j in 1:length(allModels)) { # for each species
    df <- data.frame(ModelName = NA,
                     AIC = NA,
                     diffFromBest = NA)
    for(m in 1:length(allModels[[j]])){
        # remove all models with NaN in standard error of model output
        result <- tryCatch(
            {
                summary(allModels[[j]][[m]])
            },
            error = function(e) e
        )
        if (inherits(result, "error")) {
            df[m, 1:2] <- NA
        } else {
            if (is.na(summary(allModels[[j]][[m]])$state$SE[1]) == TRUE) {
                df[m, 1] <- NA
                df[m, 2] <- NA
            } else {
                df[m, 1] <- as.character(c(allModels[[j]][[m]]@formula))
                df[m, 2] <- allModels[[j]][[m]]@AIC
            }
        }

        # if (is.na(summary(allModels[[j]][[m]])$state$SE[1]) == TRUE) {
        #     df[m, 1] <- NA
        #     df[m, 2] <- NA
        #     } else {
        #         df[m, 1] <- as.character(c(allModels[[j]][[m]]@formula))
        #         df[m, 2] <- allModels[[j]][[m]]@AIC
        #     }
      
    }
    df <- df[order(df$AIC),]
    df <- df[!is.na(df$ModelName), ]
    df$diffFromBest <- df$AIC - min(df$AIC)
    modelAICs[[j]] <- df
    
    # only the best
    ANTM <- subset(df, diffFromBest <= 2)
    topModels[[j]] <- ANTM
  }
  names(topModels) <- speciesNames
  
  masterTopModels[[i]] <- topModels
  
  # take these best models and put them into a list
  bestModsFitLists <- list()
  bestModsOutputs <- list()
  for (j in 1:length(topModels)){ # for each species
    speciesBestMods <- list()
    for(m in 1:nrow(topModels[[j]])) {
      test <- occu(formula(topModels[[j]]$ModelName[m]), ufoMasterList[[i]][[j]])
      
      speciesBestMods[[m]] <- occu(formula(topModels[[j]]$ModelName[m]), ufoMasterList[[i]][[j]], 
                                   control = 10000, 
                                   starts = c(rep(0, length(test@opt$par)))) 
    } 
    bestModsFitLists[[j]] <- fitList(speciesBestMods)
    bestModsOutputs[[j]] <- speciesBestMods
  }
  names(bestModsFitLists) <- speciesNames
  names(bestModsOutputs) <- speciesNames
  
  masterBestModsFitLists[[i]] <- bestModsFitLists
  masterBestModsOutputs[[i]] <- bestModsOutputs
  
  # just the #1 model for each species
  bestOfTheBest <- list()
  for (j in 1:length(speciesNames)) {
    bestOfTheBest[[j]] <- occu(formula(topModels[[j]]$ModelName[1]), ufoMasterList[[i]][[j]])
    
  }
  names(bestOfTheBest) <- speciesNames
  
  masterBestofTheBest[[i]] <- bestOfTheBest
  
}

# added 10/08/24 because Community isn't in best models, so can't look per community
# calculate the global model for each species (minus distToWater since that wasn't in any best model)
globalModels <- list()
for (j in 1:length(speciesNames)){ # for each species
      test <- occu(~Community + DaysEffortScaled ~RainfallScaled + Ag10KMScaled + NatArea10KMScaled + DistToComm + TemperatureScaled + Community, 
                    ufoMasterList[[1]][[j]]) # i = 1 when all communities are together

      globalModels[[j]] <- occu(~Community + DaysEffortScaled ~RainfallScaled + Ag10KMScaled + NatArea10KMScaled + DistToComm +
                    TemperatureScaled + Community, 
                    ufoMasterList[[1]][[j]], 
                                   control = 10000, 
                                   starts = c(rep(0, length(test@opt$par)))) 
      #globalModelsFitList[[j]] <- fitList(global = globalModels[[j]])
      print(paste("Done with", j, "out of", length(speciesNames), "species :)"))
}
names(globalModels) <- speciesNames
globalModelsNested <- list(Global = globalModels) # need to nest so it matches the dimensions of masterBestModsFitLists moving forwards

communityOnlyModels <- list()
for (j in 1:length(speciesNames)) { # for each species
    test <- occu(
        ~Community ~ Community,
        ufoMasterList[[1]][[j]]
    ) # i = 1 when all communities are together

    communityOnlyModels[[j]] <- occu(
        ~Community ~ Community,
        ufoMasterList[[1]][[j]],
        control = 10000,
        starts = c(rep(0, length(test@opt$par)))
    )
    # globalModelsFitList[[j]] <- fitList(global = globalModels[[j]])
    print(paste("Done with", j, "out of", length(speciesNames), "species :)"))
}
names(communityOnlyModels) <- speciesNames
communityOnlyModelsNested <- list(Global = communityOnlyModels) # need to nest so it matches the dimensions of masterBestModsFitLists moving forwards


# if proceeding with forced global models:
if (useGlobalModels == "YES") {
    masterBestModsFitLists <- list()
    masterBestModsFitLists <- globalModelsNested
    names(masterBestModsFitLists) <- communities
} 



################################################################################
################################################################################
############################# MODEL AVERAGING ##################################
################################################################################
################################################################################

# load the modelAveraging.r function I built
source("Global/modelAveraging.r")

# model average all the model outputs and get coefficients for all covariates
# make an empty list with the same length as the number of species
modelAverages <- list()

# model average the best models
for (j in 1:length(speciesNames)) {
    avg <- ManualModelAverage(masterBestModsOutputs[[1]][[j]], speciesName = speciesNames[j])
    modelAverages[[j]] <- avg
}
names(modelAverages) <- speciesNames

# save model averages to r object
save(modelAverages, file = "Global/Data/R Objects/modelAverages.RData")
load("Global/Data/R Objects/modelAverages.RData") # loads 'modelAverages'

# rbind all the $occupancy and $detection from each list
# bind_rows automatically fills missing columns with NA
occupancyList <- lapply(modelAverages, function(x) x$occupancy)
detectionList <- lapply(modelAverages, function(x) x$detection)

modelAveragesDF <- bind_rows(occupancyList)
detectionModelAverages <- bind_rows(detectionList)

# pivot the data frame to long format
# convert all non-Species columns to numeric (forcing "<NA>" strings to real NAs)
modelAveragesDF_clean <- modelAveragesDF %>%
    mutate(across(-Species, ~ as.numeric(as.character(.))))

# reshape the data
longDF <- modelAveragesDF_clean %>%
    pivot_longer(
        cols = -Species,
        names_to = "full_param",
        values_to = "value"
    ) %>%
    mutate(
        param = str_replace(full_param, "(SE|Upper|Lower)$", ""),
        type = case_when(
            str_detect(full_param, "SE$") ~ "se",
            str_detect(full_param, "Upper$") ~ "upper",
            str_detect(full_param, "Lower$") ~ "lower",
            TRUE ~ "estimate"
        )
    ) %>%
    select(Species, param, type, value) %>%
    pivot_wider(
        names_from = type,
        values_from = value
    )

# check
head(longDF)

# make significant column
longDF <- longDF %>%
    mutate(significant = ifelse(lower > 0 | upper < 0, "Significant", "Not Significant")) 
longDF$significant <- factor(longDF$significant, levels = c("Significant", "Not Significant"))


save(
    list = c(
        "masterBestModsFitLists", "speciesNames", "commonNames", "masterBestModsOutputs",
        "siteCovariate", "stations", "modelAveragesDF_clean"
    ),
    file = "Global/Data/R Objects/masterBestModsFitLists.RData"
)



################################################################################
########################## PLOTTING EFFECT SIZES ################################
################################################################################

# parameters:
# "CommunityZabalo"    "RainfallScaled"     "NatArea10KMScaled" 
# [4] "TemperatureScaled"  "DistToComm"         "CommunityRemolino" 
# [7] "CommunitySinangoe"  "CommunitySan.Pablo" "CommunitySiona"    

unique(longDF$param) # verify order
paramLabels <- data.frame(
    paramLabels = c(
        "Territory: Zábalo",
        "Avg. rainfall (kg/m²/s, scaled)",
        "% nat. area (10 km, scaled)",
        "Avg. temp. (°C, scaled)",
        "Distance to community (km)",
        "Territory: Remolino",
        "Territory: Sinangoe",
        "Territory: San Pablo",
        "Territory: Siona"
    ),
    param = c(
        "CommunityZabalo", "RainfallScaled", "NatArea10KMScaled",
        "TemperatureScaled", "DistToComm", "CommunityRemolino",
        "CommunitySinangoe", "CommunitySan.Pablo", "CommunitySiona"
    )
)

longDF <- merge(longDF, paramLabels, by = "param")



ggplot(longDF, aes(x = estimate, y = Species, color = significant)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    #geom_text(aes(x = upper + 3, label = significant), size = 8, hjust = 0, color = "red") + # nudge asterisk to the right
    facet_wrap(~paramLabels, scales = "free_x") +
    labs(x = "Effect size (95% CI)", y = NULL) +
    scale_color_manual(breaks = c("Significant", "Not Significant"), values = c("Significant" = "darkred", "Not Significant" = "gray40")) +
    theme_minimal() +
    theme( 
        text = element_text(family = "Times", colour = "black"),
        legend.title = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "italic"),
        axis.title.x = element_text(size = 14),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA)
    )

# save model averaged effect sizes 
ggsave("Global/Figures/SingleSpeciesModeling/ModelAveragedEffectSizes.png", 
width = 10, height = 8)


# extract covariate base names
cols <- colnames(modelAveragesDF)
covariates <- unique(gsub("(SE|Upper|Lower)$", "", grep("Upper|Lower|SE", cols, value = TRUE)))

# function to pivot one covariate at a time for estimates and CIs
get_covariate_df <- function(df, covariate) {
    est_col <- covariate
    lower_col <- paste0(covariate, "Lower")
    upper_col <- paste0(covariate, "Upper")

    # Use numeric vectors with NA if columns missing
    Estimate <- if (est_col %in% names(df)) as.numeric(df[[est_col]]) else rep(NA_real_, nrow(df))
    Lower <- if (lower_col %in% names(df)) as.numeric(df[[lower_col]]) else rep(NA_real_, nrow(df))
    Upper <- if (upper_col %in% names(df)) as.numeric(df[[upper_col]]) else rep(NA_real_, nrow(df))

    tibble(
        Species = df$Species,
        Covariate = covariate,
        Estimate = Estimate,
        Lower = Lower,
        Upper = Upper
    )
}

# build long table for all covariates
cov_list <- lapply(covariates, get_covariate_df, df = modelAveragesDF)
long_df <- bind_rows(cov_list)

# format confidence intervals and bold estimates where CI excludes zero
long_df <- long_df %>%
    mutate(
        # Format CI as ["lower", "upper"], NA if missing
        ConfInt = ifelse(
            !is.na(Lower) & !is.na(Upper),
            paste0('[', round(Lower, 4), ', ', round(Upper, 4), ']'),
            NA_character_
        ),
        # Check if CI excludes zero to bold Estimate
        BoldEstimate = ifelse(
            !is.na(Lower) & !is.na(Upper) & (Lower > 0 | Upper < 0),
            paste0(round(Estimate, 4)),
            as.character(round(Estimate, 4))
        )
    )

# prep species column to only show species name once in bold (first occurrence)
long_df <- long_df %>%
    arrange(Species, Covariate) %>%
    group_by(Species) %>%
    mutate(
        Species_bold = ifelse(row_number() == 1, paste0(Species[1]), "")
    ) %>%
    ungroup()


# select and rename columns for the final table
formatted_table <- long_df %>%
    mutate(
        `95% CI` = sprintf("[%.2f, %.2f]", Lower, Upper),
        Estimate = ifelse(Lower > 0 | Upper < 0,
            sprintf("<b>%.3f</b>", Estimate),
            sprintf("%.3f", Estimate)
        )
    ) %>%
    select(Species_bold, Covariate, Estimate, `95% CI`)

final_table <- formatted_table %>%
    select(
        Species = Species_bold,
        Covariate,
        Estimate = Estimate,
        ConfidenceInterval = `95% CI`
    ) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "-", .))) %>%
    mutate(across(everything(), ~ ifelse(grepl("NA", .), "-", .)))

# print with kable
kbl(final_table,
    col.names = c("Species", "Covariate", "Estimate", "95% CI"),
    escape = FALSE
) %>%
    kable_classic(full_width = TRUE, html_font = "TimesNewRoman") %>%
    column_spec(1, bold = TRUE, italic = TRUE) %>%
    save_kable(file = "Global/Figures/SingleSpeciesModeling/modelAveragesTable.png", zoom = 1.5)
















################################################################################
########################## ESTIMATE CALCULATING ################################
################################################################################

# output lists:
masterEstimatedParameters <- list() # df of occ/det and their SE/range
masterUnmarkedPredOcc <- list() # unmarked::predict() output
masterUnmarkedPredDet <- list() # unmarked::predict() output
perCovPerCommDetection <- list()
perCovPerCommOccupancy <- list()
perCovDetection <- list()
perCovOccupancy <- list()
masterGlobalDetectionEstimates <- list()
masterGlobalOccupancyEstimates <- list()

for (i in 1:length(communities)) {
    estimatedParameters <- list()
    unmarkedPredOcc <- list()
    unmarkedPredDet <- list()

    # stations info
    stations <- read.csv(paste0(communities[i], "/Data/", communitiesAbrv[i], "StationsFormatted.csv"))
    cameraRecords <- read.csv(paste0(communities[i], "/Data/", communitiesAbrv[i], "IndependentRecordsFormatted.csv"))
    # replace all Mazama species with Mazama sp.
    cameraRecords$Species <- gsub("Mazama americana", "Mazama sp.", cameraRecords$Species)
    cameraRecords$Species <- gsub("Mazama nemorivaga", "Mazama sp.", cameraRecords$Species)
    cameraRecords$Species <- gsub("Mazama gouazoubira", "Mazama sp.", cameraRecords$Species) # replace all Mazama species with Mazama sp.

    # import site covariates for each community
    if (communities[i] == "Sinangoe") {
        siteCovariate <- data.frame(DistToComm = scale(stations$Community / 1000)) # site covariates (scaled)
        df <- data.frame(
            DistToComm = siteCovariate$DistToComm
        )

        for (j in 1:length(speciesNames)) { # for each species
            # state = occupancy
            unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                type = "state", new = df, appendData = TRUE
            )
            # det = detection
            bestModel <- occu(
                formula(masterTopModels[[i]][[j]]$ModelName[1]),
                ufoMasterList[[i]][[j]]
            )
            # nullDetPred <- backTransform(bestModel, "det")

            # null
            nullModel <- occu(~1 ~ 1, ufoMasterList[[1]][[j]])
            nullDetPred <- backTransform(nullModel, "det")
            nullOccuPred <- backTransform(nullModel, "state")

            # average
            estimatedParameters[[j]] <- data.frame(
                avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                nullOccupancy = nullOccuPred@estimate,
                nullOccupancySE = SE(nullOccuPred),
                nullDetection = nullDetPred@estimate,
                nullDetectionSE = SE(nullDetPred),
                detectionRange = NA
            )
        }
    } else if (communities[i] == "Zabalo") {
        load("Zabalo/Data/R Objects/siteCovs2018.RData") # loads 'siteCovariate'
        df <- data.frame(
            Habitat = siteCovariate$Habitat,
            HuntingIntensity = siteCovariate$HuntingIntensity,
            Trail.Distance = siteCovariate$Trail.Distance,
            Effort = siteCovariate$Effort
        )

        for (j in 1:length(speciesNames)) { # for each species
            # state = occupancy
            unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                type = "state", new = df, appendData = TRUE
            )
            # det = detection
            unmarkedPredDet[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                type = "det", new = df, appendData = TRUE
            )

            # average
            estimatedParameters[[j]] <- data.frame(
                avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                avgDetection = mean(unmarkedPredDet[[j]]$Predicted),
                avgDetectionSE = mean(unmarkedPredDet[[j]]$SE),
                detectionRange = range(unmarkedPredDet[[j]]$Predicted)
            )
        }
    } else if (communities[i] == "Siona" | communities[i] == "Siekopai") {
        # for Siona and Siekopai who don't have covariates
        siteCovariate <- NULL # no covariates for remaining communities
        df <- NULL # not sure what to do with this since we don't have covariates for other communities

        # calculate occupancy and detection estimates
        for (j in 1:length(speciesNames)) {
            bestModel <- occu(~1 ~ 1, ufoMasterList[[i]][[j]])
            nullOccPred <- backTransform(bestModel, "state")
            nullDetPred <- backTransform(bestModel, "det")

            # average
            estimatedParameters[[j]] <- data.frame(
                avgOccupancy = nullOccPred@estimate,
                avgOccupancySE = SE(nullOccPred),
                occupancyRange = NA,
                avgDetection = nullDetPred@estimate,
                avgDetectionSE = SE(nullDetPred),
                detectionRange = NA
            )
        }
    } else if (communities[i] == "Global") {
        covariates <- c(
            "Community", "NatArea10KMScaled", "RainfallScaled",
            "DistToWater", "TemperatureScaled", "DistToComm"
        )
        siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
        siteCovariate$Rainfall <- siteCovariate$Rainfall * 1000
        allCommunities <- unique(siteCovariate$Community)
        N <- 50
        dfTemplate <- data.frame(
            Community = rep(allCommunities, each = N), 
            RainfallScaled = mean(siteCovariate$RainfallScaled),
            NatArea10KMScaled = mean(siteCovariate$NatArea10KMScaled),
            DistToWater = mean(siteCovariate$DistToWater),
            TemperatureScaled = mean(siteCovariate$TemperatureScaled),
            DistToComm = mean(siteCovariate$DistToComm),
            DaysEffortScaled = mean(siteCovariate$DaysEffortScaled)
        )

        for (j in 1:length(speciesNames)) { # for each species
            dfEdited <- dfTemplate
            ########## average overall
            # state = occupancy
            unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                type = "state", new = dfTemplate, appendData = TRUE
            )
            # det = detection
            unmarkedPredDet[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                type = "det", new = dfTemplate, appendData = TRUE
            )

            # null
            nullModel <- occu(~1 ~ 1, ufoMasterList[[i]][[j]])
            nullDetPred <- backTransform(nullModel, "det")
            nullOccuPred <- backTransform(nullModel, "state")

            # average
            estimatedParameters[[j]] <- data.frame(
                avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                avgDetection = mean(unmarkedPredDet[[j]]$Predicted),
                avgDetectionSE = mean(unmarkedPredDet[[j]]$SE),
                nullOccupancy = nullOccuPred@estimate,
                nullOccupancySE = SE(nullOccuPred),
                nullDetection = nullDetPred@estimate,
                nullDetectionSE = SE(nullDetPred),
                detectionRange = NA
            )

            # for (m in 1:length(allCommunities)) {
            #     dfEdited <- dfTemplate
            #     dfEdited[, "Community"] <- allCommunities[m]

            ######## prediction per covariate
            for (k in 1:length(covariates)) {
                dfEdited <- dfTemplate
                # dfEdited[, "Community"] <- allCommunities[m]
                covariateInQuestion <- covariates[k]

                # prediction DF
                if (covariateInQuestion == "Community") {
                    dfEdited[, covariateInQuestion] <- rep(rep(unique(siteCovariate$Community),
                        each = N / length(unique(siteCovariate$Community))
                    ), times = length(allCommunities))
                } else {
                    # dfEdited[, "Community"] <- allCommunities[m]
                    dfEdited[, covariateInQuestion] <- seq(min(siteCovariate[, covariateInQuestion]),
                        max(siteCovariate[, covariateInQuestion]),
                        length.out = N
                    )
                }

                # state = occupancy
                df <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                    type = "state", new = dfEdited, appendData = TRUE
                )
                df$PredictedCovariate <- covariateInQuestion
                # df$CommunityHeldConstant <- allCommunities[m]
                df$Species <- speciesNames[j]
                perCovOccupancy[[k]] <- df

                # det = detection
                df <- unmarked::predict(masterBestModsFitLists[[i]][[j]],
                    type = "det", new = dfEdited, appendData = TRUE
                )
                df$PredictedCovariate <- covariateInQuestion
                # df$CommunityHeldConstant <- allCommunities[m]
                df$Species <- speciesNames[j]
                perCovDetection[[k]] <- df
            }
            names(perCovOccupancy) <- covariates
            names(perCovDetection) <- covariates

            perCovPerCommOccupancy <- perCovOccupancy
            perCovPerCommDetection <- perCovDetection

            # }

            # names(perCovPerCommOccupancy) <- allCommunities
            # names(perCovPerCommDetection) <- allCommunities

            masterGlobalOccupancyEstimates[[j]] <- perCovPerCommOccupancy
            masterGlobalDetectionEstimates[[j]] <- perCovPerCommDetection
        }
        names(masterGlobalOccupancyEstimates) <- speciesNames
        names(masterGlobalDetectionEstimates) <- speciesNames
    }

    names(estimatedParameters) <- speciesNames

    masterEstimatedParameters[[i]] <- estimatedParameters
    masterUnmarkedPredOcc[[i]] <- unmarkedPredOcc
    masterUnmarkedPredDet[[i]] <- unmarkedPredDet
}

names(masterEstimatedParameters) <- communities
names(masterUnmarkedPredOcc) <- communities
names(masterUnmarkedPredDet) <- communities




################################################################################
############################## PLOT ESTIMATES ##################################
################################################################################

communitiesAccent <- gsub(pattern = "Zabalo", replacement = "Zábalo", communities)
speciesNames
#masterEstimatedParameters[[1]]

# make a dataframe for facilitated plotting
estimates <- data.frame(Community = rep(communitiesAccent, each = length(speciesNames)),
                        Species = rep(speciesNames, times = length(communitiesAccent)),
                        CommonNames = rep(commonNames, times = length(communitiesAccent)),
                        nullOccupancy = NA,
                        nullOccupancySE = NA,
                        nullDetection = NA,
                        nullDetectionSE = NA,
                        avgOccupancy = NA,
                        avgOccupancySE = NA,
                        avgOccupancySEManual = NA,
                        avgOccupancySD = NA,
                        avgDetection = NA,
                        avgDetectionSE = NA,
                        avgDetectionSEManual = NA,
                        avgDetectionSD = NA)
for (i in 1:length(communitiesAccent)) {
  for (j in 1:length(speciesNames)){
    row <- which((estimates$Community == communitiesAccent[i]) & (estimates$Species == speciesNames[j]))
    estimates[row, "avgOccupancy"] <- masterEstimatedParameters[[i]][[j]]$avgOccupancy[1]
    estimates[row, "avgOccupancySE"] <- masterEstimatedParameters[[i]][[j]]$avgOccupancySE[1]
    estimates[row, "nullOccupancy"] <- masterEstimatedParameters[[i]][[j]]$nullOccupancy[1]
    estimates[row, "nullOccupancySE"] <- masterEstimatedParameters[[i]][[j]]$nullOccupancySE[1]
    estimates[row, "nullDetection"] <- masterEstimatedParameters[[i]][[j]]$nullDetection[1]
    estimates[row, "nullDetectionSE"] <- masterEstimatedParameters[[i]][[j]]$nullDetectionSE[1]
    #estimates[row, "avgOccupancySEManual"] <- sd(masterUnmarkedPredOcc[[i]][[j]]$Predicted)/sqrt(nrow(masterUnmarkedPredOcc[[i]][[j]]))
    #estimates[row, "avgOccupancySD"] <- sd(masterUnmarkedPredOcc[[i]][[j]]$Predicted)
    estimates[row, "avgDetection"] <- masterEstimatedParameters[[i]][[j]]$avgDetection[1]
    estimates[row, "avgDetectionSE"] <- masterEstimatedParameters[[i]][[j]]$avgDetectionSE[1]
    #estimates[row, "avgDetectionSEManual"] <- sd(masterUnmarkedPredDet[[i]][[j]]$Predicted)/sqrt(nrow(masterUnmarkedPredDet[[i]][[j]]))
    #estimates[row, "avgDetectionSD"] <- sd(masterUnmarkedPredDet[[i]][[j]]$Predicted)
  }
}
estimates$Species <- factor(estimates$Species, levels = speciesNames) # so plotting doesn't alphabetize species



################################################################################
################################################################################
######################## REQUIRES MANUAL INPUT!!! ##############################
################################################################################
################################################################################


# DO YOU WANT TO SAVE PLOTS WHEN RUNNING?????? THIS WILL OVERWRITE OLD PLOTS!!!!
savePlots <- "YES" # "YES" or "NO"


################################################################################
################################################################################
################################################################################
################################################################################




# FROM HERE: PLOTS WILL ONLY WORK IF SPECIES == PECCARY, BROCKET, PACA, TRUMPETER
commonNames




################################################################################
############################## PLOT ESTIMATES ##################################
################################################################################


# make a dataframe that has predictions across each covariate to facilitate plotting
covariatesMinusCommunity <- covariates[covariates != "Community"] # everything but community

# making the data frame to plot from
plottingDF <- as.data.frame(do.call(rbind, do.call(rbind, masterGlobalOccupancyEstimates)))
plottingDF <- plottingDF[plottingDF$PredictedCovariate != "Community",]
plottingDF$Community <- gsub("Zabalo", "Zábalo", x = plottingDF$Community)
plottingDF$Species <- factor(plottingDF$Species, levels = speciesNames) # so plotting doesn't alphabetize species

# order the communities by percent of natural cover
orderedCommunities <- c(siteCovariate %>% 
                          group_by(Community) %>% 
                          summarize(perc = mean(natArea10KM)) %>% 
                          arrange(desc(perc)) %>% 
                          select(Community))$Community
orderedCommunities <- gsub("Zabalo", "Zábalo", x = orderedCommunities)
plottingDF <- plottingDF %>% dplyr::filter(Community %in% orderedCommunities)
plottingDF$Community <- factor(plottingDF$Community, levels = orderedCommunities)
plottingDF$Species <- factor(plottingDF$Species, levels = speciesNames) # so plotting doesn't alphabetize species


# add species common names
crossList <- data.frame(commonNames = commonNames,
                        Species = speciesNames)
plottingDF <- merge(plottingDF, crossList, by = "Species")

# make labels per covariate
colors <- c("Zábalo" = "darkgreen", "Remolino" = "forestgreen", 
            "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3")

# make axis titles per species
peccary <- ~ atop(paste(bold("Collared peccary")), paste("(", italic("Pecari tajacu"), ")"))
brocket <- ~ atop(paste(bold("Brocket")), paste("(", italic("Mazama sp."), ")"))
paca <- ~ atop(paste(bold("Lowland paca")), paste("(", italic("Cuniculus paca"), ")"))
trumpeter <- ~ atop(paste(bold("Grey-winged trumpeter")), paste("(", italic("Psophia crepitans"), ")"))
fourEyed <- ~ atop(paste(bold("Brown four-eyed opossum")), paste("(", italic("Metachirus nudicaudatus"), ")"))
agouti <- ~ atop(paste(bold("Black agouti")), paste("(", italic("Dasyprocta fuliginosa"), ")"))
armadillo <- ~ atop(paste(bold("Nine-banded armadillo")), paste("(", italic("Dasypus novemcinctus"), ")"))
tinamou <- ~ atop(paste(bold("Great tinamou")), paste("(", italic("Tinamus major"), ")"))
opossum <- ~ atop(paste(bold("Common opossum")), paste("(", italic("Didelphis marsupialis"), ")"))
ocelot <- ~ atop(paste(bold("Ocelot")), paste("(", italic("Leopardus pardalis"), ")"))

# # rphylopic per species
# peccPic <- get_uuid(name = "Pecari tajacu", n = 1)
# brockPic <- get_phylopic(uuid = get_uuid(name = "Mazama americana", n = 1))
# pacaPic <- get_phylopic(uuid = get_uuid(name = "Cuniculus paca", n = 1))
# trumpPic <- get_phylopic(uuid = get_uuid(name = "Psophia crepitans", n = 1))
# #fourEyedPic <- get_phylopic(uuid = get_uuid(name = "Metachirus nudicaudatus", n = 1))
# agoutiPic <- get_phylopic(uuid = get_uuid(name = "Dasyprocta", n = 1))
# armadilloPic <- get_phylopic(uuid = get_uuid(name = "Dasypus novemcinctus", n = 1))
# tinamouPic <- get_phylopic(uuid = get_uuid(name = "Tinamus major", n = 1))
# #opossumPic <- get_phylopic(uuid = get_uuid(name = "Didelphis", n = 1))
# ocelotPic <- get_phylopic(uuid = get_uuid(name = "Leopardus pardalis", n = 1))

# plot null occupancy
# ensure order matches speciesNames vector
estimates$Species <- factor(estimates$Species, levels = speciesNames)

# URLs
peccaryURL <- "https://images.phylopic.org/images/44fb7d4f-6d59-432b-9583-a87490259789/raster/1024x610.png?v=183ff0ad631"
brocketURL <- "https://images.phylopic.org/images/b5f40112-0cb8-4994-aa70-28ac97ccb83f/raster/901x1024.png?v=186bb13accc"
pacaURL <- "https://images.phylopic.org/images/414b0720-a160-4bce-b060-2eb9675fc1c8/raster/1024x559.png?v=17d5f8f529f"
trumpeterURL <- "https://images.phylopic.org/images/feb8e7c3-483d-4f78-8d9e-5618e96102e7/raster/565x1024.png?v=1985dfa2256"
agoutiURL <- "https://images.phylopic.org/images/30fe5e82-8127-4cbb-9c3f-c64a379376a8/raster/1024x642.png?v=178d877cf4e"
armadilloURL <- "https://images.phylopic.org/images/5d59b5ce-c1dd-40f6-b295-8d2629b9775e/raster/1024x493.png?v=13571ec0b6a"
tinamouURL <- "https://images.phylopic.org/images/446debec-ca63-4882-801f-beaf479887d5/raster/1024x773.png?v=1464c44f011"
ocelotURL <- "https://images.phylopic.org/images/2fc7bbbf-8ca7-48fb-8495-d351cd3b1f99/raster/1024x409.png?v=176942952dd"

# URLs and desired heights
sil_urls <- c(
    peccaryURL, brocketURL, pacaURL, trumpeterURL,
    agoutiURL, armadilloURL, tinamouURL, ocelotURL
)
sil_heights <- c(0.10, 0.125, 0.10, 0.13, 0.10, 0.10, 0.125, 0.10)

# fade helper: multiply existing alpha by `alpha` and write to a temp PNG
fade_images <- function(urls, alpha = 0.2) {
    vapply(urls, function(u) {
        img <- image_read(u)
        img <- image_fx(img, expression = paste0("u*", alpha), channel = "alpha")
        f <- tempfile(fileext = ".png")
        image_write(img, path = f, format = "png")
        f
    }, FUN.VALUE = character(1))
}

faded_pngs <- fade_images(sil_urls, alpha = 0.2)

sil_df <- data.frame(
    Species = speciesNames,
    y = 0.05,
    image = faded_pngs,
    height = sil_heights
)

p <- ggplot(estimates, aes(
    x = Species,
    y = nullOccupancy
)) +
    geom_point(size = 1.5) +
    geom_errorbar(
        aes(
            ymin = nullOccupancy - nullOccupancySE,
            ymax = nullOccupancy + nullOccupancySE
        ),
        width = 0.15, linewidth = .5
    ) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "black") +
    scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter, agouti, armadillo, tinamou, ocelot)) +
    labs(x = "Species", y = "Null occupancy probability (SE)") +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.60),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)
    ) +
    geom_image(
        data = sil_df,
        aes(x = Species, y = y, image = image, size = height),
        alpha = 0.7,
        by = "height" ) +
    scale_size_identity()

p

# save it
if (savePlots == "YES") {
  ggsave(filename = paste0(communities, "/Figures/SingleSpeciesModeling/", 
                           communitiesAbrv, "nullOccupancyEstimates.png"), 
         width = 8, height = 4)
  }

############ PLOT BY COMMUNITY ##############
# making the data frame to plot from
communityPlottingDF <- as.data.frame(do.call(rbind, do.call(rbind, masterGlobalOccupancyEstimates)))
communityPlottingDF <- communityPlottingDF[communityPlottingDF$PredictedCovariate == "Community",]
communityPlottingDF$Community <- gsub("Zabalo", "Zábalo", x = communityPlottingDF$Community)
communityPlottingDF$Species <- factor(communityPlottingDF$Species, levels = speciesNames) # so plotting doesn't alphabetize species
communityPlottingDF$Community <- factor(communityPlottingDF$Community,
    levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)
communityPlottingDF$Species <- factor(communityPlottingDF$Species, levels = speciesNames) # so plotting doesn't alphabetize species

communityPlottingDF <- communityPlottingDF %>% distinct()



### Check for significance
## See what is statistically significant
model <- lm(Predicted ~ Species, data = communityPlottingDF)
summary(model)
test <- aov(Predicted ~ Species, data = communityPlottingDF)
summary(test)
TukeyHSD(test)



# plot it
dodge <- position_dodge(width = 0.3)
p <- ggplot(communityPlottingDF, aes(
    x = Species,
    y = Predicted,
    color = Community
)) +
    geom_point(aes(color = Community), position = dodge, size = 1.5) +
    geom_errorbar(
        aes(
            ymin = lower,
            ymax = upper,
            color = Community
        ),
        position = dodge, width = 0.15, linewidth = .5
    ) +
    # scale_color_manual(values = c("darkorange", "royalblue", "green3", "yellow3")) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter, agouti, armadillo, tinamou, ocelot)) +
    labs(x = "Species", y = "Predicted occupancy probability (95% CI)") +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.60),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        #legend.position = "none",
        panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)
    ) +
    geom_image(
        data = sil_df,
        aes(x = Species, y = y, image = image, size = height),
        alpha = 0.7,
        inherit.aes = FALSE, # prevents inheriting color=Community
        by = "height" ) +
    scale_size_identity()

# plot with the animal silhouettes
p

# save it
if (savePlots == "YES") {
    ggsave(
        filename = paste0(
            communities, "/Figures/SingleSpeciesModeling/",
            communitiesAbrv, "OccupancyEstimates.png"
        ),
        width = 8, height = 4
    )
}



################################################################################
################################################################################
######################## REQUIRES MANUAL INPUT!!! ##############################
################################################################################
################################################################################


# DO YOU WANT TO SAVE PLOTS WHEN RUNNING?????? THIS WILL OVERWRITE OLD PLOTS!!!!
savePlots <- "YES" # "YES" or "NO"


################################################################################
################################################################################
################################################################################
################################################################################









################################################################################
############################# PLOT PREDICTIONS #################################
################################################################################

# make a plot faceted by species for each covariate
# with prediction lines divided by Community if it was a covariate in the model
for (i in 1:length(covariatesMinusCommunity)){
  # covariate in question
  cov <- covariatesMinusCommunity[i]
  df <- plottingDF[plottingDF$PredictedCovariate == cov,]
  covariateInQuestion <- df[,cov]
  
  # x axis labels for each covariate
  if (covariatesMinusCommunity[i] == "RainfallScaled") {
    xlabel <- expression(paste("Rainfall (g*", m^{-2}, s^{-1}, ", scaled)"))
  } else if (covariatesMinusCommunity[i] == "TemperatureScaled"){
    xlabel <- "Temperature (°C, scaled)"
  } else if (covariatesMinusCommunity[i] == "DistToWater") {
    xlabel <- "Distance to water source (m)"
  } else if (covariatesMinusCommunity[i] == "NatArea10KMScaled"){
    xlabel <- "Percent natural area within 10 km (scaled)"
  } else if (covariatesMinusCommunity[i] == "DistToComm") {
    xlabel <- "Distance to a community (km)"
  } else if (covariatesMinusCommunity[i] == "DaysEffortScaled") {
    xlabel <- "Camera trap effort (days, scaled)"
  } else if (covariatesMinusCommunity[i] == "Ag10KMScaled") {
    xlabel <- "Percent agricultural area within 10 km (scaled)"
  }
  
  # xlab:
  # rainfall: expression(paste("Rainfall (g*", m^{-2}, s^{-1}, ", scaled)"))
  # percentNatural: Percent natural area within 25 km
  # Temperature: Temperature (°C, scaled)
  # DistToWater: Distance to water (m)
  # DistToComm: Distance to a community (km)
  
  # actually plot it
  ggplot(df, aes(x = covariateInQuestion, 
                 y = Predicted, color = Community)) +
    geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), 
                alpha = 0.2, color = NA) +
    geom_line(aes(x = covariateInQuestion, y = Predicted)) +
    ylab(expression(paste("Occupancy probability estimate (", psi, ")"))) +
    facet_wrap(~commonNames, nrow = 3, ncol = 4) +
    xlab(xlabel) +
    coord_cartesian(ylim = c(0,1), 
                    xlim = c(min(covariateInQuestion),
                             max(covariateInQuestion))) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(text = element_text(family = "Times", colour = "black"),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          #axis.title.x = element_blank(), 
          legend.position="top"
    )
  #theme(plot.title = element_text(hjust = 0.5))
  
  # save it
  if (savePlots == "YES") {
    ggsave(filename = paste0(communities, "/Figures/SingleSpeciesModeling/", 
                             cov, "Prediction.png"), 
           width = 8, height = 4)
  }
  
}
nDaysGroupedPerSpecies





# ################################################################################
# ########################### PLOT MAP OF STATIONS ###############################
# ################################################################################

# # load data
# data("world")
# stations <- read.csv("Global/Data/AllStationsFormatted.csv")
# stations <- stations %>%
#   select(c(Station, gps_x, gps_y, CommunityName)) %>%
#   distinct()
# stations$CommunityName <- gsub("Zabalo", "Zábalo", x = stations$CommunityName)

# # order the communities by percent of natural cover
# stations <- stations %>% dplyr::filter(CommunityName %in% orderedCommunities)
# stations$Community <- factor(stations$CommunityName, levels=orderedCommunities)
# stations$CommunityName <- NULL

# # only highlight Ecuador
# SA <- c("ecuador", "bolivia", "brazil", "chile", "colombia", "argentina", "guyana", "paraguay", "peru", "suriname", "uruguay", "venezuela")
# mapColors <- rep("white", length(SA))
# mapColors[2] <- "lightyellow"
# colors <- c(
#     "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
#     "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
# )

# # option A with data("world")
# worldEdited <- world
# worldEdited$Ecu <- ifelse(worldEdited$name_long == "Ecuador", "A", "B")

# # option B with 'maps' package
# mappy <- maps::map("worldHires", SA)
# mappy_sf <- st_as_sf(mappy, crs = st_crs(worldEdited), fill = FALSE)
# mappy_sf$Ecu <- ifelse(mappy_sf$ID == "Ecuador", "A", "B")

# # option A: less detail in Ecuador border
# ecuadorMap <- worldEdited %>%
#   dplyr::filter(continent == "South America") %>%
#   ggplot() +
#   geom_sf(aes(fill=Ecu)) +
#   geom_point(data = stations, aes(gps_x, gps_y, fill = Community), pch = 21) +
#   #scale_fill_manual(values = c("lightyellow", "white")) +
#   scale_fill_manual(name = "Community",
#                     values = c(colors, "A" = "lightyellow", "B" = "white"), 
#                     breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")) +
#   coord_sf(default_crs = sf::st_crs(4326), xlim = c(-150, -37)) + 
#   #guides(fill = "none") +
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         text = element_text(family = "Times", colour = "black"),
#         axis.text = element_text(colour = "black"))

# # option B: more detail in Ecuador border but can't get Ecuador to highlight yellow?
# ecuadorMap <- mappy_sf %>%
#   ggplot() +
#   geom_sf(aes(fill=Ecu), lwd = 0.5) +
#   geom_point(data = stations, aes(gps_x, gps_y, fill = Community), pch = 21) +
#   #scale_fill_manual(values = c("lightyellow", "white")) +
#   scale_fill_manual(name = "Community",
#                     values = c(colors, "A" = "lightyellow", "B" = "white"), 
#                     breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona"),
#                     guide = guide_legend(override.aes = list(shape = 21, size = 6.5, fill = colors) )) +
#   coord_sf(default_crs = sf::st_crs(4326), xlim = c(-150, -37)) + 
#   #guides(fill = "none") +
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         text = element_text(family = "Times", colour = "black"),
#         axis.text = element_text(colour = "black"))

# # set bounding box to magnify
# coords <- as.matrix(stations[,c("gps_x","gps_y")])
# e <- as.vector(ext(coords)) 
# e["xmin"] <- e["xmin"] - 0.1
# e["ymin"] <- e["ymin"] - 0.5
# e["xmax"] <- e["xmax"] + 0.1
# e["ymax"] <- e["ymax"] + 0.5

# # plot map inlay
# together <- ecuadorMap + geom_magnify(from = e, 
#                           to = c(xmin = -150, xmax = -85, ymin = -45, ymax = 3))
# together

# if (savePlots == "YES") {
#   ggsave(plot = together, filename = paste0(communities, "/Figures/mapInlayWithSites.png"), 
#          width = 7, height = 5)
# }




################################################################################
######################### PLOT PERCENT NATURAL AREA ############################
################################################################################

# # plot percent natural area within 25 km with SD
# siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
# natStats <- siteCovariate %>% 
#   group_by(Community) %>%
#   summarize(avgNat = mean(percentNatural), sdNat = sd(percentNatural))
# natStats$Community <- gsub("Zabalo", "Zábalo", x = natStats$Community)

# # order the communities by percent of natural cover
# natStats <- natStats %>% dplyr::filter(Community %in% orderedCommunities)
# natStats$Community <- factor(natStats$Community, levels=orderedCommunities)

# # plot it
# ggplot(natStats, aes(x = Community, y = avgNat, fill = Community)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin = avgNat - sdNat, ymax = avgNat + sdNat), width = 0.2) +
#   ylab("Percent natural area within 25 km") +
#   scale_fill_manual(values = colors) +
#   ylim(c(0,1)) +
#   theme_bw()+
#   theme(text = element_text(family = "Times", colour = "black"),
#         axis.text = element_text(colour = "black"))
# # save it
# if (savePlots == "YES") {
#   ggsave(filename = paste0(communities, "/Figures/percentNatArea25km.png"), 
#          width = 7, height = 5)
# }

# # zoomed
# ggplot(natStats, aes(x = Community, y = avgNat, fill = Community)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin = avgNat - sdNat, ymax = avgNat + sdNat), width = 0.2) +
#   ylab("Percent natural area within 25 km") +
#   scale_fill_manual(values = colors) +
#   coord_cartesian(ylim = c(0.825,1)) +
#   theme_bw() +
#   theme(text = element_text(family = "Times", colour = "black"),
#         axis.text = element_text(colour = "black"))
# # save it
# if (savePlots == "YES") {
#   ggsave(filename = paste0(communities, "/Figures/percentNatArea25kmZoomed.png"), 
#          width = 7, height = 5)
# }






################################################################################
####################### SITE ABUNDANCE AND DIVERSITY ###########################
########################### AND CAMERA TRAP INFO ###############################
################################################################################

# load data
Data <- read.csv("Global/Data/AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("Global/Data/AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
# load site covariates for all communities
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
siteCovariate$Rainfall <- siteCovariate$Rainfall*1000 # convert to grams/m^2/s  


# replace all Mazama species with Mazama sp.


# remove all unknown species
noUnknownsSGE <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Sinangoe"),]
noUnknownsZAB <- Data[(Data$Species != "NAN NAN" & Data$Species != "NA NA" & Data$CommunityName == "Zabalo"),]
noUnknownsSNA <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Siona"),]
noUnknownsSPA <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "San Pablo"),]
noUnknownsREM <- Data[(Data$Species != "N/D N/D" & Data$CommunityName == "Remolino"),]

# camera trap info
cameraInfo <- Data %>% 
  filter(Species != "N/D N/D" & Species != "NAN NAN" & Species != "NA NA") %>%
  group_by(CommunityName) %>%
  mutate(StartDate = min(DateTimeOriginal),
         EndDate = max(as.Date(DateTimeOriginal)),
         Year = year(DateTimeOriginal)) %>%
  summarise(OperatingDays = round(as.numeric(max(DateTimeOriginal)-
                                               min(DateTimeOriginal))),
            StartDate = min(DateTimeOriginal),
            EndDate = max(DateTimeOriginal),
            numberOfStations = length(unique(Station))) 
# numberOfCamerasPerStation = round(length(unique(CameraName))/length(unique(Station)))

# format the dates into a legible format
cameraInfo$StartDate <- format(cameraInfo$StartDate, "%Y-%m-%d")
cameraInfo$EndDate <- format(cameraInfo$EndDate, "%Y-%m-%d")

# put cameraInfo so that it is in the following order, Zabalo, Remolino, Sinangoe, San Pablo, then Siona
cameraInfo <- cameraInfo[c(5,1,3,2,4),] # order by natural area
cameraInfo$CommunityName <- gsub(cameraInfo$CommunityName, pattern = "Zabalo", replacement = "Zábalo")

# per station
# Sinangoe
siteDiversitySGE <- noUnknownsSGE %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySGE <- siteDiversitySGE %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Zábalo
siteDiversityZAB <- noUnknownsZAB %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversityZAB <- siteDiversityZAB %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Siona
siteDiversitySNA <- noUnknownsSNA %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySNA <- siteDiversitySNA %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# San Pablo
siteDiversitySPA <- noUnknownsSPA %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversitySPA <- siteDiversitySPA %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# Remolino
siteDiversityREM <- noUnknownsREM %>%
  group_by(Station, Species) %>%
  summarise(abundance = n()) 

siteDiversityREM <- siteDiversityREM %>%
  group_by(Station) %>%
  summarise(N=sum(abundance),
            shannonDiversity = -sum((abundance/sum(abundance))*log(abundance/sum(abundance))),
            simpsonDiversity = 1-sum((abundance/sum(abundance))^2))

# per community
wholeDiversitySGE <- noUnknownsSGE %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Sinangoe", 
         NatArea10KM = mean(siteCovariate$NatArea10KM[siteCovariate$Community == "Sinangoe"]),
         OperatingDays = round(as.numeric(max(noUnknownsSGE$DateTimeOriginal)-
                                            min(noUnknownsSGE$DateTimeOriginal))))
wholeDiversityZAB <- noUnknownsZAB %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Zábalo", 
         NatArea10KM = round(mean(siteCovariate$NatArea10KM[siteCovariate$Community == "Zabalo"]), 3), 
         OperatingDays = round(as.numeric(max(noUnknownsZAB$DateTimeOriginal)-
                                            min(noUnknownsZAB$DateTimeOriginal))))
wholeDiversitySNA <- noUnknownsSNA %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Siona", 
         NatArea10KM = round(mean(siteCovariate$NatArea10KM[siteCovariate$Community == "Siona"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsSNA$DateTimeOriginal)-
                                            min(noUnknownsSNA$DateTimeOriginal))))
wholeDiversitySPA <- noUnknownsSPA %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "San Pablo", 
         NatArea10KM = round(mean(siteCovariate$NatArea10KM[siteCovariate$Community == "San Pablo"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsSPA$DateTimeOriginal)-
                                            min(noUnknownsSPA$DateTimeOriginal))))

wholeDiversityREM <- noUnknownsREM %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Remolino", 
         NatArea10KM = round(mean(siteCovariate$NatArea10KM[siteCovariate$Community == "Remolino"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsREM$DateTimeOriginal)-
                                            min(noUnknownsREM$DateTimeOriginal))))

# abundance and diversity for all communities
communityAbundance <- rbind(wholeDiversityZAB, wholeDiversitySPA, wholeDiversityREM, 
                            wholeDiversitySGE, wholeDiversitySNA)

communityDiversity <- communityAbundance %>%
  group_by(Community, NatArea10KM) %>%
  summarise(nIndiv=sum(abundance),
            nSpecies = length(unique(Species)),
            OperatingDays = mean(OperatingDays),
            shannonIndex = round(-sum((abundance/sum(abundance))*log(abundance/sum(abundance))), 3),
            simpsonIndex = round(1-sum((abundance/sum(abundance))^2), 3)) 
communityDiversity$NatArea10KM <- round(communityDiversity$NatArea10KM, 3)
communityDiversity <- arrange(communityDiversity, desc(NatArea10KM))
communityDiversity$Community <- factor(communityDiversity$Community, 
                                       levels = communityDiversity$Community)
communityDiversity
write.csv(communityDiversity, "Global/Data/CommunityDiversityAbundance.csv")



######### TABLE OF THE NUMBER OF DETECTIONS PER SPECIES IN DATA
nDetectionsPerSpecies <- Data %>%
  filter(Species != "N/D N/D" & Species != "NAN NAN" & Species != "NA NA") %>%
  group_by(Species) %>%
  summarise(nDetections = n()) %>%
  arrange(desc(nDetections))
nDetectionsPerSpecies %>% 
    kbl(col.names = c("Species", "Number of Detections")) %>%
    kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
    row_spec(c(1, 2, 3, 5, 6, 7, 8, 13,18,21, 32), bold = TRUE) %>%
    kableExtra::save_kable(file = "Global/Figures/numberOfDetectionsPerSpecies.png", zoom = 10)










################################################################################
################################################################################
######################## REQUIRES MANUAL INPUT!!! ##############################
################################################################################
################################################################################


# DO YOU WANT TO SAVE PLOTS WHEN RUNNING?????? THIS WILL OVERWRITE OLD PLOTS!!!!
savePlots <- "YES" # "YES" or "NO"


################################################################################
################################################################################
################################################################################
################################################################################








################################################################################
################################ MAKE TABLES ###################################
################################################################################

# peek
head(communityDiversity)
head(cameraInfo)

# save it
if (savePlots == "YES") {
  # table with diversity and abundance information
  kbl(communityDiversity, col.names = c("Territory", "Proportion of Natural Area (within 10 km)", 
                                        "Number of Detections", "Number of Species",
                                        "Number of Sampling Days",
                                        "Shannon Diversity Index", "Simpson Diversity Index")) %>%
    kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
    kableExtra::save_kable(file = "Global/Figures/communityDiversityAbundance.png", zoom = 1.5)
  
  # table with just diversity information
  kbl(communityDiversity[,c("Community", "NatArea10KM", "OperatingDays", 
                            "shannonIndex", "simpsonIndex")], 
      col.names = c("Territory", "Proportion of Natural Area (within 10 km)", "Number of Sampling Days", 
                    "Shannon Diversity Index", "Simpson Diversity Index")) %>%
    kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
    save_kable(file = "Global/Figures/communityDiversitySummary.png", zoom = 2)
  
  # table with camera trap information
  kbl(cameraInfo[,1:4], col.names = c("Territory", "Number of Sampling Days", 
                                "Sampling Start Date", "Sampling End Date")) %>%
    kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
    save_kable(file = "Global/Figures/siteInfo.png",
               zoom = 10)

}

# Make a table with average occupancy and detection estimates from above
if (communities == "Global" & savePlots == "YES"){
    # table with occupancy and detection estimates

    # merge the estimates with nDaysGroupedPerSpecies to show the number of days grouped
    estimatesWithDays <- merge(estimates, nDaysGroupedPerSpecies, by = "CommonNames")

    estimatesCut <- estimatesWithDays |> 
      dplyr::select(CommonNames, DaysGrouped, avgOccupancy, avgOccupancySE, avgDetection, avgDetectionSE) |>
      arrange(desc(avgOccupancy)) |>
      dplyr::mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCut, col.names = c("Species", "Days per Detection Occasion",
                                 "Average Occupancy Estimate", "Occupancy SE", 
                                 "Average Detection Estimate", "Detection SE")) %>%
      kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
      kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/occupancyDetectionEstimates.png", zoom = 10)
    
    estimatesCutNull <- estimatesWithDays |>
        dplyr::select(CommonNames, DaysGrouped, nullOccupancy, nullOccupancySE, nullDetection, nullDetectionSE) |>
        arrange(desc(nullOccupancy)) |>
        dplyr::mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCutNull, col.names = c(
        "Species", "Days per Detection Occasion",
        "Null Occupancy Estimate", "Null Occupancy SE",
        "Null Detection Estimate", "Null Detection SE"
    )) %>%
        kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
        kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/nullOccupancyDetectionEstimates.png", zoom = 10)

    estimatesCutBoth <- estimatesWithDays |>
        dplyr::select(CommonNames, DaysGrouped, avgOccupancy, avgOccupancySE, nullOccupancy, nullOccupancySE,
                      avgDetection, avgDetectionSE, nullDetection, nullDetectionSE) |>
        arrange(desc(avgOccupancy)) |>
        dplyr::mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCutBoth, col.names = c(
        "Species", "Days per Detection Occasion",
        "Average Modeled Occupancy Estimate", "Average Modeled Occupancy SE",
        "Null Occupancy Estimate", "Null Occupancy SE",
        "Average Modeled Detection Estimate", "Average Modeled Detection SE",
        "Null Detection Estimate", "Null Detection SE"
    )) %>%
        kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
        kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/occupancyDetectionEstimates_Both.png", zoom = 10)

}


# make a table with the covariates included in the top models for each species
if (communities == "Global" & savePlots == "YES" & all(speciesNames == names(masterTopModels[[1]])) & length(speciesNames) == 8){
    
    for (i in 1:length(speciesNames)){
        # make a column with the species and how many best models
        masterTopModels[[1]][[i]]$Species <- names(masterTopModels[[1]])[i]
        masterTopModels[[1]][[i]]$nModels <- nrow(masterTopModels[[1]][[i]])
    }

    # combine all the species into one dataframe
    bestModelsDF <- do.call(rbind.data.frame, masterTopModels[[1]])
    rownames(bestModelsDF) <- NULL

    # format said dataframe
    bestModelsDF$AIC <- round(bestModelsDF$AIC, 3)
    bestModelsDF$diffFromBest <- round(bestModelsDF$diffFromBest, 3)

    # https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html#Group_rows_via_labeling
    nModelsPerSpecies <- c(bestModelsDF %>%
        distinct(Species, nModels) %>%
        select(nModels))$nModels

    kbl(bestModelsDF[, c("ModelName", "AIC", "diffFromBest")],
        col.names = c("Model", "AIC", "Delta AIC"),
        # caption = names(masterTopModels[[1]])[i],
    ) %>% # AIC weight?
        kable_classic(full_width = TRUE, html_font = "TimesNewRoman") %>%
        pack_rows(speciesNames[1],
            start_row = 1,
            end_row = nModelsPerSpecies[1]
        ) %>%
        pack_rows(speciesNames[2],
            start_row = nModelsPerSpecies[1] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2]
        ) %>%
        pack_rows(speciesNames[3],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3]
        ) %>%
        pack_rows(speciesNames[4],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4]
        ) %>%
        pack_rows(speciesNames[5],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5]
        ) %>%
        pack_rows(speciesNames[6],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + nModelsPerSpecies[6]
        ) %>%
        pack_rows(speciesNames[7],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + nModelsPerSpecies[6] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + nModelsPerSpecies[6] + nModelsPerSpecies[7]
        ) %>%
        pack_rows(speciesNames[8],
            start_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + nModelsPerSpecies[6] + nModelsPerSpecies[7] + 1,
            end_row = nModelsPerSpecies[1] + nModelsPerSpecies[2] + nModelsPerSpecies[3] + nModelsPerSpecies[4] + nModelsPerSpecies[5] + nModelsPerSpecies[6] + nModelsPerSpecies[7] + nModelsPerSpecies[8]
        ) %>%
        kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/AllBestModelsTable.png", zoom = 2)


}




################################################################################
######################## YEAR AS A COVARIATE ANALYSIS ##########################
################################################################################
# Goal: Test whether sampling year (2018 vs. 2022) influences occupancy or
# detection probability, to justify treating both periods as a single season

yearModels <- list()
yearAICResults <- list()

for (j in 1:length(speciesNames)) {
    ufo <- ufoMasterList[[1]][[j]] # [[1]] = Global (the only community here)

    # cast Year as a factor so it is treated as categorical (2018 vs. 2022)
    ufo@siteCovs$Year <- as.factor(ufo@siteCovs$Year)

    # four candidate models: null, year on detection only,
    # year on occupancy only, year on both
    yearModelFormulas <- list(
        "Null"               = "~ 1 ~ 1",
        "Year (detection)"   = "~ Year ~ 1",
        "Year (occupancy)"   = "~ 1 ~ Year",
        "Year (both)"        = "~ Year ~ Year"
    )

    speciesYearMods <- list()
    for (m in 1:length(yearModelFormulas)) {
        test <- occu(formula(yearModelFormulas[[m]]), ufo)
        speciesYearMods[[m]] <- occu(
            formula(yearModelFormulas[[m]]), ufo,
            control = 10000,
            starts = rep(0, length(test@opt$par))
        )
    }
    names(speciesYearMods) <- names(yearModelFormulas)
    yearModels[[j]] <- speciesYearMods

    # AIC table for this species
    aicDF <- data.frame(
        Species = commonNames[j],
        Model = names(yearModelFormulas),
        AIC = round(sapply(speciesYearMods, function(m) m@AIC), 3),
        stringsAsFactors = FALSE
    )
    aicDF$DeltaAIC <- round(aicDF$AIC - min(aicDF$AIC), 3)
    yearAICResults[[j]] <- aicDF

    print(paste("Year model done for", commonNames[j]))
}
names(yearModels) <- speciesNames

# combine all species into dataframe
yearAICTable <- do.call(rbind, yearAICResults)
rownames(yearAICTable) <- NULL

# sort by species and then AIC
yearAICTable <- yearAICTable %>%
    arrange(Species, AIC)

# save it
if (savePlots == "YES") {
    n <- length(yearModelFormulas) # number of models per species (always 4)
    cumN <- cumsum(rep(n, length(commonNames)))

    kbl(yearAICTable[, c("Model", "AIC", "DeltaAIC")],
        col.names = c("Model", "AIC", "\u0394AIC")
    ) %>%
        kable_classic(full_width = TRUE, html_font = "TimesNewRoman") %>%
        pack_rows(commonNames[1], 1, cumN[1]) %>%
        pack_rows(commonNames[2], cumN[1] + 1, cumN[2]) %>%
        pack_rows(commonNames[3], cumN[2] + 1, cumN[3]) %>%
        pack_rows(commonNames[4], cumN[3] + 1, cumN[4]) %>%
        pack_rows(commonNames[5], cumN[4] + 1, cumN[5]) %>%
        pack_rows(commonNames[6], cumN[5] + 1, cumN[6]) %>%
        pack_rows(commonNames[7], cumN[6] + 1, cumN[7]) %>%
        pack_rows(commonNames[8], cumN[7] + 1, cumN[8]) %>%
        kableExtra::save_kable(
            file = "Global/Figures/SingleSpeciesModeling/yearModelAICComparison.png",
            zoom = 2
        )
}


################################################################################
#################### FULL MODEL SELECTION WITH YEAR ###########################
################################################################################
# Year added as a candidate covariate on occ

# Year as a factor in every species' UFO (Global community = index 1)
ufoListYear <- ufoMasterList[[1]]
for (j in 1:length(ufoListYear)) {
    ufoListYear[[j]]@siteCovs$Year <- as.factor(ufoListYear[[j]]@siteCovs$Year)
}

# candidate covariate sets; same as original Global section, plus Year
match_detection_year <- c("Community", "DaysEffortScaled", "Year")
match_occupancy_year <- c(
    "Community", "RainfallScaled", "NatArea10KMScaled",
    "DistToWater", "TemperatureScaled", "DistToComm", "Year"
)


################################################################################
###################### STAGE 1: BEST DETECTION MODEL ##########################
################################################################################

# every combination of detection covariates
combos <- sapply(seq(length(match_detection_year)), function(k) {
    as.list(as.data.frame(combn(x = match_detection_year, m = k)))
})
combos <- unlist(combos, recursive = FALSE)

forms <- sapply(combos, function(x) paste("~", paste(x, collapse = "+")))
detectionFormulas_year <- as.vector(c(forms, "~ 1"))

tempDF <- data.frame(detection = detectionFormulas_year, occupancy = "~ 1")
allDetectionFormulas_year <- paste(tempDF$detection, tempDF$occupancy, sep = " ")

bestDetectionModels_year <- data.frame(species = commonNames, bestDetectionModel = NA)

for (j in 1:length(speciesNames)) {
    temp <- data.frame(Model = allDetectionFormulas_year, AIC = NA)

    for (m in 1:length(allDetectionFormulas_year)) {
        test <- occu(formula(allDetectionFormulas_year[[m]]), ufoListYear[[j]])
        mod <- occu(formula(allDetectionFormulas_year[[m]]), ufoListYear[[j]],
            control = 10000,
            starts  = rep(0, length(test@opt$par))
        )
        temp$AIC[m] <- mod@AIC
    }

    bestDetectionModels_year$bestDetectionModel[j] <-
        detectionFormulas_year[which(temp$AIC == min(temp$AIC))]

    # override: keep ocelot on null detection
    if (commonNames[j] == "Ocelot") {
        bestDetectionModels_year$bestDetectionModel[j] <- "~ 1"
    }

    print(paste("Year detection stage done for", commonNames[j]))
}


################################################################################
###################### STAGE 2: BEST OCCUPANCY MODEL ##########################
################################################################################

# every combination of occupancy covariates
combos <- sapply(seq(length(match_occupancy_year)), function(k) {
    as.list(as.data.frame(combn(x = match_occupancy_year, m = k)))
})
combos <- unlist(combos, recursive = FALSE)

forms <- sapply(combos, function(x) paste("~", paste(x, collapse = "+")))
occupancyFormulas_year <- as.vector(c(forms, "~ 1"))

# pair best detection formula with every occupancy formula
occupancyModelsList_year <- list()
for (j in 1:length(speciesNames)) {
    df <- data.frame(
        detection = bestDetectionModels_year$bestDetectionModel[j],
        occupancy = occupancyFormulas_year
    )
    occupancyModelsList_year[[j]] <- paste(df$detection, df$occupancy, sep = " ")
}

# all occupancy models
allModels_year <- list()
for (j in 1:length(speciesNames)) {
    # override: keep ocelot on null
    if (commonNames[j] == "Ocelot") {
        occupancyModelsList_year[[j]] <- "~ 1 ~ 1"
    }

    occupancyMods <- list()
    for (m in 1:length(occupancyModelsList_year[[j]])) {
        test <- occu(formula(occupancyModelsList_year[[j]][m]), ufoListYear[[j]])
        occupancyMods[[m]] <- occu(formula(occupancyModelsList_year[[j]][m]), ufoListYear[[j]],
            control = 10000,
            starts  = rep(0, length(test@opt$par))
        )
    }
    names(occupancyMods) <- 1:length(occupancyMods)
    allModels_year[[j]] <- occupancyMods

    print(paste("Year occupancy stage done for", j, "out of", length(speciesNames), "species :)"))
}


################################################################################
########################## TOP MODEL SELECTION #################################
################################################################################

modelAICs_year <- list()
topModels_year <- list()

for (j in 1:length(allModels_year)) {
    df <- data.frame(ModelName = NA, AIC = NA, diffFromBest = NA)

    for (m in 1:length(allModels_year[[j]])) {
        result <- tryCatch(summary(allModels_year[[j]][[m]]), error = function(e) e)

        if (inherits(result, "error")) {
            df[m, 1:2] <- NA
        } else {
            if (is.na(summary(allModels_year[[j]][[m]])$state$SE[1])) {
                df[m, 1:2] <- NA
            } else {
                df[m, 1] <- as.character(c(allModels_year[[j]][[m]]@formula))
                df[m, 2] <- allModels_year[[j]][[m]]@AIC
            }
        }
    }

    df <- df[order(df$AIC), ]
    df <- df[!is.na(df$ModelName), ]
    df$diffFromBest <- df$AIC - min(df$AIC)
    modelAICs_year[[j]] <- df

    ANTM <- subset(df, diffFromBest <= 2)
    topModels_year[[j]] <- ANTM
}
names(topModels_year) <- speciesNames


################################################################################
####################### DID YEAR MAKE IT INTO TOP MODELS? #####################
################################################################################

yearInTopModels <- data.frame(
    Species = commonNames,
    YearInTopModel = sapply(topModels_year, function(df) {
        any(grepl("Year", df$ModelName))
    }),
    ModelsWithYear = sapply(topModels_year, function(df) {
        matched <- df$ModelName[grepl("Year", df$ModelName)]
        if (length(matched) == 0) "None" else paste(matched, collapse = "; ")
    }),
    stringsAsFactors = FALSE
)

print(yearInTopModels)


################################################################################
################################## TABLE #######################################
################################################################################

if (savePlots == "YES") {
    # full top-model AIC table (same style as AllBestModelsTable)
    for (j in 1:length(speciesNames)) {
        topModels_year[[j]]$Species <- commonNames[j]
        topModels_year[[j]]$nModels <- nrow(topModels_year[[j]])
    }

    bestModelsDF_year <- do.call(rbind.data.frame, topModels_year)
    rownames(bestModelsDF_year) <- NULL
    bestModelsDF_year$AIC <- round(bestModelsDF_year$AIC, 3)
    bestModelsDF_year$diffFromBest <- round(bestModelsDF_year$diffFromBest, 3)

    nModelsPerSpecies_year <- c(bestModelsDF_year %>%
        distinct(Species, nModels) %>%
        select(nModels))$nModels
    cumN <- cumsum(nModelsPerSpecies_year)

    kbl(bestModelsDF_year[, c("ModelName", "AIC", "diffFromBest")],
        col.names = c("Model", "AIC", "\u0394AIC")
    ) %>%
        kable_classic(full_width = TRUE, html_font = "TimesNewRoman") %>%
        pack_rows(speciesNames[1], 1, cumN[1]) %>%
        pack_rows(speciesNames[2], cumN[1] + 1, cumN[2]) %>%
        pack_rows(speciesNames[3], cumN[2] + 1, cumN[3]) %>%
        pack_rows(speciesNames[4], cumN[3] + 1, cumN[4]) %>%
        pack_rows(speciesNames[5], cumN[4] + 1, cumN[5]) %>%
        pack_rows(speciesNames[6], cumN[5] + 1, cumN[6]) %>%
        pack_rows(speciesNames[7], cumN[6] + 1, cumN[7]) %>%
        pack_rows(speciesNames[8], cumN[7] + 1, cumN[8]) %>%
        kableExtra::save_kable(
            file = "Global/Figures/SingleSpeciesModeling/AllBestModelsTable_withYear.png",
            zoom = 2
        )

    # compact summary table: which species had Year in a top model
    kbl(yearInTopModels,
        col.names = c("Species", "Year in Top Model?", "Model(s) Containing Year")
    ) %>%
        kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
        kableExtra::save_kable(
            file = "Global/Figures/SingleSpeciesModeling/yearInTopModelsSummary.png",
            zoom = 2
        )
}


# chi-squared test: correlation between community and sampling year
year_community_table <- table(siteCovariate$Community, siteCovariate$Year)
print(year_community_table)
chisq.test(year_community_table)

################################################################################
################################# PERMANOVA ####################################
################################################################################
# GOAL: test whether animal community composition at camera stations differs
# between sampling years (zabalo 2018 vs. all other communities 2022) 

require(vegan)

ufoGlobal <- ufoMasterList[[1]]

detRateMatrix <- sapply(seq_along(ufoGlobal), function(j) {
    y <- ufoGlobal[[j]]@y
    rowSums(y, na.rm = TRUE) / rowSums(!is.na(y))
})
colnames(detRateMatrix) <- commonNames
rownames(detRateMatrix) <- rownames(ufoGlobal[[1]]@y)

# pull year and community from the site covariates already attached to the ufos
siteMeta <- ufoGlobal[[1]]@siteCovs[, c("Community", "Year")]
siteMeta$Year <- as.factor(siteMeta$Year)
siteMeta$Community <- as.factor(siteMeta$Community)

# remove stations with NA or all-zero detection rates
badRows <- which(rowSums(is.na(detRateMatrix)) > 0 | rowSums(detRateMatrix) == 0)
detRateMatrix <- detRateMatrix[-badRows, ]
siteMeta <- siteMeta[-badRows, ]

# confirm row alignment
stopifnot(nrow(detRateMatrix) == nrow(siteMeta))

# model 1: year alone
set.seed(123)
permanova_year <- adonis2(detRateMatrix ~ Year,
    data = siteMeta,
    method = "bray",
    permutations = 999
)
print(permanova_year)

# model 2: year after conditioning on community
# asks whether year explains residual variation beyond community identity
set.seed(123)
permanova_year_conditioned <- adonis2(detRateMatrix ~ Community + Year,
    data = siteMeta,
    method = "bray",
    permutations = 999
)
print(permanova_year_conditioned)

# homogeneity of dispersions check (assumption of adonis2)
# if groups differ in spread rather than location, results need caveating
bray_dist <- vegdist(detRateMatrix, method = "bray")
dispersion_year <- betadisper(bray_dist, siteMeta$Year)
permutest(dispersion_year, permutations = 999)

# nmds
set.seed(704)
nmds <- metaMDS(detRateMatrix, distance = "bray", k = 2, trymax = 100)
cat("NMDS stress:", nmds$stress, "\n")
# stress < 0.10 = good; < 0.20 = acceptable for ecological data

# build a data frame for ggplot
nmdsScores <- as.data.frame(scores(nmds, display = "sites"))
nmdsScores$Year <- siteMeta$Year
nmdsScores$Community <- siteMeta$Community

# compute group centroids per year
yearCentroids <- nmdsScores %>%
    group_by(Year) %>%
    summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2))

nmdsPlot <- ggplot(nmdsScores, aes(x = NMDS1, y = NMDS2, colour = Year, fill = Year)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, linetype = "dashed") +
    geom_point(size = 3, alpha = 0.8) +
    geom_point(
        data = yearCentroids, aes(x = NMDS1, y = NMDS2),
        shape = 23, size = 5, colour = "black", stroke = 1
    ) +
    scale_colour_manual(
        values = c("2018" = "#E69F00", "2022" = "#56B4E9"),
        name = "Sampling year"
    ) +
    scale_fill_manual(
        values = c("2018" = "#E69F00", "2022" = "#56B4E9"),
        name = "Sampling year"
    ) +
    annotate("text",
        x = min(nmdsScores$NMDS1),
        y = max(nmdsScores$NMDS2),
        label = paste0("Stress = ", round(nmds$stress, 3)),
        hjust = 0, size = 4
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12)
    )

print(nmdsPlot)

if (savePlots == "YES") {
    ggsave("Global/Figures/SingleSpeciesModeling/nmdsByYear.png",
        plot = nmdsPlot,
        width = 6,
        height = 5,
        dpi = 300
    )
}

# format the adonis2 output into a supplementary table
if (savePlots == "YES") {
    permanovaDF <- data.frame(
        Term = c(
            "Year", "Residual", "Total",
            "Community", "Year (conditioned)", "Residual", "Total"
        ),
        Model = c(rep("Year only", 3), rep("Community + Year", 4)),
        Df = c(permanova_year$Df, permanova_year_conditioned$Df),
        SumOfSqs = round(c(permanova_year$SumOfSqs, permanova_year_conditioned$SumOfSqs), 4),
        R2 = round(c(permanova_year$R2, permanova_year_conditioned$R2), 4),
        F = round(c(permanova_year$F, permanova_year_conditioned$F), 3),
        p.value = c(permanova_year$`Pr(>F)`, permanova_year_conditioned$`Pr(>F)`)
    )

    # replace NA in non-test rows with em-dashes
    permanovaDF[is.na(permanovaDF)] <- "\u2014"

    kbl(permanovaDF[, c("Term", "Df", "SumOfSqs", "R2", "F", "p.value")],
        col.names = c("Term", "df", "Sum of squares", "R\u00B2", "F", "p"),
        row.names = FALSE
    ) %>%
        kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
        pack_rows("Model 1: Year only", 1, 3) %>%
        pack_rows("Model 2: Community + Year", 4, 7) %>%
        kableExtra::save_kable(
            file = "Global/Figures/SingleSpeciesModeling/permanovaResultsTable.png",
            zoom = 2
        )
}










# TIME!
#toc() 


