
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
require(rphylopic)
require(ggpubr)
require(sf)
require(spData)
require(ggmagnify)
require(terra)
require(mapdata)
require(kableExtra)
require(knitr)
require(lubridate)
require(tictoc)
require(MuMIn)
require(tidyverse)

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
    match_occupancy <- c("Community", "percentNatural", "RainfallScaled",
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
      test <- occu(~Community + DaysEffortScaled ~RainfallScaled + percentNatural + DistToComm + TemperatureScaled + Community, 
                    ufoMasterList[[1]][[j]]) # i = 1 when all communities are together
      
      globalModels[[j]] <- occu(~Community + DaysEffortScaled ~RainfallScaled + percentNatural + DistToComm +
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

  

#################################################################################
################################################################################
################################################################################
# Model averaging

# model average all the model outputs and get coefficients for all covariates
# make an empty list with the same length as the number of species
modelAverages <- list()
# model average the best models
for (j in 1:length(speciesNames)) {
    if (length(masterBestModsOutputs[[1]][[j]]) == 1) {
        modelAverages[[j]] <- masterBestModsOutputs[[1]][[j]][[1]]@estimates
    } else {
        # model average the best models
        bestModsOutputs <- masterBestModsOutputs[[1]][[j]]
        # model average the best models
        avg <- model.avg(bestModsOutputs)
        summ <- summary(avg)

        # save full model averages to list
        modelAverages[[j]] <- summ # less biased
    }
}
names(modelAverages) <- speciesNames


# save model averages to r object
save(modelAverages, file = "Global/Data/R Objects/modelAverages.RData")
load("Global/Data/R Objects/modelAverages.RData") # loads 'modelAverages'

# make a dataframe with the model averaged coefficients with each covariate as a column
modelAveragesDF <- data.frame(
    Species = names(modelAverages),
    percentNatural = NA,
    percentNaturalSE = NA,
    percentNaturalUpper = NA,
    percentNaturalLower = NA,
    RainfallScaled = NA,
    RainfallScaledSE = NA,
    RainfallScaledUpper = NA,
    RainfallScaledLower = NA,
    DistToWater = NA,
    DistToWaterSE = NA,
    DistToWaterUpper = NA,
    DistToWaterLower = NA,
    TemperatureScaled = NA,
    TemperatureScaledSE = NA,
    TemperatureScaledUpper = NA,
    TemperatureScaledLower = NA,
    DistToComm = NA,
    DistToCommSE = NA,
    DistToCommUpper = NA,
    DistToCommLower = NA,
    CommunityZabalo = NA,
    CommunityZabaloSE = NA,
    CommunityZabaloUpper = NA,
    CommunityZabaloLower = NA,
    CommunitySinangoe = NA,
    CommunitySinangoeSE = NA,
    CommunitySinangoeUpper = NA,
    CommunitySinangoeLower = NA,
    CommunitySanPablo = NA,
    CommunitySanPabloSE = NA,
    CommunitySanPabloUpper = NA,
    CommunitySanPabloLower = NA,
    CommunityRemolino = NA,
    CommunityRemolinoSE = NA,
    CommunityRemolinoUpper = NA,
    CommunityRemolinoLower = NA,
    CommunitySiona = NA,
    CommunitySionaSE = NA,
    CommunitySionaUpper = NA,
    CommunitySionaLower = NA
)
detectionModelAverages <- data.frame(
    Species = names(modelAverages),
    DaysEffortScaled = NA,
    DaysEffortScaledSE = NA,
    DaysEffortScaledUpper = NA,
    DaysEffortScaledLower = NA,
    CommunityZabalo = NA,
    CommunityZabaloSE = NA,
    CommunityZabaloUpper = NA,
    CommunityZabaloLower = NA,
    CommunitySinangoe = NA,
    CommunitySinangoeSE = NA,
    CommunitySinangoeUpper = NA,
    CommunitySinangoeLower = NA,
    CommunitySanPablo = NA,
    CommunitySanPabloSE = NA,
    CommunitySanPabloUpper = NA,
    CommunitySanPabloLower = NA,
    CommunityRemolino = NA,
    CommunityRemolinoSE = NA,
    CommunityRemolinoUpper = NA,
    CommunityRemolinoLower = NA,
    CommunitySiona = NA,
    CommunitySionaSE = NA,
    CommunitySionaUpper = NA,
    CommunitySionaLower = NA
)
for (i in 1:length(modelAverages)) {
#coefTable(summ, full = TRUE)
#confint(summ, full = TRUE)
#modelAverages[[i]]$coefficients["full", ]

    if (names(modelAverages)[i] != "Leopardus pardalis") {
        coefficients <- as.data.frame(gsub(" ", "", rbind(
            names(modelAverages[[i]]$coefficients["full", ]),
            modelAverages[[i]]$coefficients["full", ])
        ))
        coefficients <- rbind(
            coefficients,
            as.data.frame(confint(modelAverages[[i]], full = TRUE))[, "2.5 %"], # lower
            as.data.frame(confint(modelAverages[[i]], full = TRUE))[, "97.5 %"] # upper
        ) 
        
        for (j in 1:ncol(coefficients)){
            # if the column is a psi column
            if (grepl("psi", coefficients[1, j]) == TRUE) {
                # if Community is a covariate in ANY column, then put the value from "psi(Int)" for Zabalo
                if (any(grepl("Community", coefficients[1, j]))) {
                    modelAveragesDF[i, "CommunityZabalo"] <- coefficients[2, "psi(Int)"]
                    modelAveragesDF[i, "CommunityZabaloLower"] <- coefficients[3, "psi(Int)"]
                    modelAveragesDF[i, "CommunityZabaloUpper"] <- coefficients[4, "psi(Int)"]
                    modelAveragesDF[i, "CommunitySinangoe"] <- coefficients[2, "psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySinangoeLower"] <- coefficients[3, "psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySinangoeUpper"] <- coefficients[4, "psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySanPablo"] <- coefficients[2, "psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunitySanPabloLower"] <- coefficients[3, "psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunitySanPabloUpper"] <- coefficients[4, "psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunityRemolino"] <- coefficients[2, "psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunityRemolinoLower"] <- coefficients[3, "psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunityRemolinoUpper"] <- coefficients[4, "psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunitySiona"] <- coefficients[2, "psi(CommunitySiona)"]
                    modelAveragesDF[i, "CommunitySionaLower"] <- coefficients[3, "psi(CommunitySiona)"]
                    modelAveragesDF[i, "CommunitySionaUpper"] <- coefficients[4, "psi(CommunitySiona)"]
                } else if (any(grepl("Rainfall", coefficients[1, j]))) {
                    modelAveragesDF[i, "RainfallScaled"] <- coefficients[2, "psi(RainfallScaled)"]
                    modelAveragesDF[i, "RainfallScaledLower"] <- coefficients[3, "psi(RainfallScaled)"]
                    modelAveragesDF[i, "RainfallScaledUpper"] <- coefficients[4, "psi(RainfallScaled)"]
                } else if (any(grepl("percentNatural", coefficients[1, j]))) {
                    modelAveragesDF[i, "percentNatural"] <- coefficients[2, "psi(percentNatural)"]
                    modelAveragesDF[i, "percentNaturalLower"] <- coefficients[3, "psi(percentNatural)"]
                    modelAveragesDF[i, "percentNaturalUpper"] <- coefficients[4, "psi(percentNatural)"]
                } else if (any(grepl("DistToWater", coefficients[1, j]))) {
                    modelAveragesDF[i, "DistToWater"] <- coefficients[2, "psi(DistToWater)"]
                    modelAveragesDF[i, "DistToWaterLower"] <- coefficients[3, "psi(DistToWater)"]
                    modelAveragesDF[i, "DistToWaterUpper"] <- coefficients[4, "psi(DistToWater)"]
                } else if (any(grepl("TemperatureScaled", coefficients[1, j]))) {
                    modelAveragesDF[i, "TemperatureScaled"] <- coefficients[2, "psi(TemperatureScaled)"]
                    modelAveragesDF[i, "TemperatureScaledLower"] <- coefficients[3, "psi(TemperatureScaled)"]
                    modelAveragesDF[i, "TemperatureScaledUpper"] <- coefficients[4, "psi(TemperatureScaled)"]
                } else if (any(grepl("DistToComm", coefficients[1, j]))) {
                    modelAveragesDF[i, "DistToComm"] <- coefficients[2, "psi(DistToComm)"]
                    modelAveragesDF[i, "DistToCommLower"] <- coefficients[3, "psi(DistToComm)"]
                    modelAveragesDF[i, "DistToCommUpper"] <- coefficients[4, "psi(DistToComm)"]
                } else next

            } else if (grepl("p", coefficients[1, j]) == TRUE) {

                # if Community is a covariate, then put the value from "psi(Int)" for Zabalo
                if (any(grepl("Community", coefficients[1, j]))) {
                    detectionModelAverages[i, "CommunityZabalo"] <- coefficients[2, "p(Int)"]
                    detectionModelAverages[i, "CommunityZabaloLower"] <- coefficients[3, "p(Int)"]
                    detectionModelAverages[i, "CommunityZabaloUpper"] <- coefficients[4, "p(Int)"]
                    detectionModelAverages[i, "CommunitySinangoe"] <- coefficients[2, "p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySinangoeLower"] <- coefficients[3, "p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySinangoeUpper"] <- coefficients[4, "p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySanPablo"] <- coefficients[2, "p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunitySanPabloLower"] <- coefficients[3, "p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunitySanPabloUpper"] <- coefficients[4, "p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunityRemolino"] <- coefficients[2, "p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunityRemolinoLower"] <- coefficients[3, "p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunityRemolinoUpper"] <- coefficients[4, "p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunitySiona"] <- coefficients[2, "p(CommunitySiona)"]
                    detectionModelAverages[i, "CommunitySionaLower"] <- coefficients[3, "p(CommunitySiona)"]
                    detectionModelAverages[i, "CommunitySionaUpper"] <- coefficients[4, "p(CommunitySiona)"]
                } else if (any(grepl("DaysEffortScaled", coefficients[1, j]))) {
                    detectionModelAverages[i, "DaysEffortScaled"] <- coefficients[2, "p(DaysEffortScaled)"]
                    detectionModelAverages[i, "DaysEffortScaledLower"] <- coefficients[3, "p(DaysEffortScaled)"]
                    detectionModelAverages[i, "DaysEffortScaledUpper"] <- coefficients[4, "p(DaysEffortScaled)"]
                } else next

            } 
        }


    } else  if (names(modelAverages)[i] == "Leopardus pardalis"){ # if it is the ocelot
        model <- modelAverages[[i]]
        coefficients <- c(coef(model@estimates$state), coef(model@estimates$det))
        lowerState <- as.data.frame(confint(model@estimates$state))[, "0.025"]
        names(lowerState) <- rownames(as.data.frame(confint(model@estimates$state)))
        lowerDet <- as.data.frame(confint(model@estimates$det))[, "0.025"]
        names(lowerDet) <- rownames(as.data.frame(confint(model@estimates$det)))
        lowers <- c(lowerState, lowerDet)
        upperState <- as.data.frame(confint(model@estimates$state))[, "0.975"]
        names(upperState) <- rownames(as.data.frame(confint(model@estimates$state)))
        upperDet <- as.data.frame(confint(model@estimates$det))[, "0.975"]
        names(upperDet) <- rownames(as.data.frame(confint(model@estimates$det)))
        uppers <- c(upperState, upperDet)
        SEs <- c(SE(model@estimates$state), SE(model@estimates$det))
        
        # for every coefficient
        for (j in 1:length(coefficients)) {
            # if the column is a psi column
            if (grepl("psi", names(coefficients)[j]) == TRUE) {
                # if Community is a covariate in ANY column, then put the value from "psi(Int)" for Zabalo
                if (any(grepl("Community", names(coefficients)[j]))) {
                    modelAveragesDF[i, "CommunityZabalo"] <- coefficients["psi(Int)"]
                    modelAveragesDF[i, "CommunityZabaloUpper"] <- uppers["psi(Int)"]
                    modelAveragesDF[i, "CommunityZabaloLower"] <- lowers["psi(Int)"]
                    modelAveragesDF[i, "CommunitySinangoe"] <- coefficients["psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySinangoeUpper"] <- uppers["psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySinangoeLower"] <- lowers["psi(CommunitySinangoe)"]
                    modelAveragesDF[i, "CommunitySanPablo"] <- coefficients["psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunitySanPabloUpper"] <- uppers["psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunitySanPabloLower"] <- lowers["psi(CommunitySan Pablo)"]
                    modelAveragesDF[i, "CommunityRemolino"] <- coefficients["psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunityRemolinoUpper"] <- uppers["psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunityRemolinoLower"] <- lowers["psi(CommunityRemolino)"]
                    modelAveragesDF[i, "CommunitySiona"] <- coefficients["psi(CommunitySiona)"]
                    modelAveragesDF[i, "CommunitySionaUpper"] <- uppers["psi(CommunitySiona)"]
                    modelAveragesDF[i, "CommunitySionaLower"] <- lowers["psi(CommunitySiona)"]
                } else if (any(grepl("RainfallScaled", names(coefficients)[j]))) {
                    modelAveragesDF[i, "RainfallScaled"] <- coefficients["psi(RainfallScaled)"]
                    modelAveragesDF[i, "RainfallScaledUpper"] <- uppers["psi(RainfallScaled)"]
                    modelAveragesDF[i, "RainfallScaledLower"] <- lowers["psi(RainfallScaled)"]
                } else if (any(grepl("percentNatural", names(coefficients)[j]))) {
                    modelAveragesDF[i, "percentNatural"] <- coefficients["psi(percentNatural)"]
                    modelAveragesDF[i, "percentNaturalUpper"] <- uppers["psi(percentNatural)"]
                    modelAveragesDF[i, "percentNaturalLower"] <- lowers["psi(percentNatural)"]
                } else if (any(grepl("DistToWater", names(coefficients)[j]))) {
                    modelAveragesDF[i, "DistToWater"] <- coefficients["psi(DistToWater)"]
                    modelAveragesDF[i, "DistToWaterUpper"] <- uppers["psi(DistToWater)"]
                    modelAveragesDF[i, "DistToWaterLower"] <- lowers["psi(DistToWater)"]
                } else if (any(grepl("TemperatureScaled", names(coefficients)[j]))) {
                    modelAveragesDF[i, "TemperatureScaled"] <- coefficients["psi(TemperatureScaled)"]
                    modelAveragesDF[i, "TemperatureScaledUpper"] <- uppers["psi(TemperatureScaled)"]
                    modelAveragesDF[i, "TemperatureScaledLower"] <- lowers["psi(TemperatureScaled)"]
                } else if (any(grepl("DistToComm", names(coefficients)[j]))) {
                    modelAveragesDF[i, "DistToComm"] <- coefficients["psi(DistToComm)"]
                    modelAveragesDF[i, "DistToCommUpper"] <- uppers["psi(DistToComm)"]
                    modelAveragesDF[i, "DistToCommLower"] <- lowers["psi(DistToComm)"]
                } else next

            } else if (grepl("p", names(coefficients)[j]) == TRUE) {

                # if Community is a covariate, then put the value from "psi(Int)" for Zabalo
                if (any(grepl("Community", names(coefficients)[j]))) {
                    detectionModelAverages[i, "CommunityZabalo"] <- coefficients["p(Int)"]
                    detectionModelAverages[i, "CommunityZabaloUpper"] <- uppers["p(Int)"]
                    detectionModelAverages[i, "CommunityZabaloLower"] <- lowers["p(Int)"]
                    detectionModelAverages[i, "CommunityZabaloSE"] <- SEs["p(Int)"]
                    detectionModelAverages[i, "CommunitySinangoe"] <- coefficients["p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySinangoeUpper"] <- uppers["p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySinangoeLower"] <- lowers["p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySinangoeSE"] <- SEs["p(CommunitySinangoe)"]
                    detectionModelAverages[i, "CommunitySanPablo"] <- coefficients["p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunitySanPabloUpper"] <- uppers["p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunitySanPabloLower"] <- lowers["p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunitySanPabloSE"] <- SEs["p(CommunitySan Pablo)"]
                    detectionModelAverages[i, "CommunityRemolino"] <- coefficients["p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunityRemolinoUpper"] <- uppers["p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunityRemolinoLower"] <- lowers["p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunityRemolinoSE"] <- SEs["p(CommunityRemolino)"]
                    detectionModelAverages[i, "CommunitySiona"] <- coefficients["p(CommunitySiona)"]
                    detectionModelAverages[i, "CommunitySionaUpper"] <- uppers["p(CommunitySiona)"]
                    detectionModelAverages[i, "CommunitySionaLower"] <- lowers["p(CommunitySiona)"]
                    detectionModelAverages[i, "CommunitySionaSE"] <- SEs["p(CommunitySiona)"]
                } else if (any(grepl("DaysEffortScaled", names(coefficients)[j]))) {
                    detectionModelAverages[i, "DaysEffortScaled"] <- coefficients["p(DaysEffortScaled)"]
                    detectionModelAverages[i, "DaysEffortScaledUpper"] <- uppers["p(DaysEffortScaled)"]
                    detectionModelAverages[i, "DaysEffortScaledLower"] <- lowers["p(DaysEffortScaled)"]
                    detectionModelAverages[i, "DaysEffortScaledSE"] <- SEs["p(DaysEffortScaled)"]
                } else next

    }

}

}
}

head(modelAveragesDF)

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



################################################################################
########################## PLOTTING EFFECT SIZES ################################
################################################################################


ggplot(longDF, aes(x = estimate, y = Species, color = significant)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    #geom_text(aes(x = upper + 3, label = significant), size = 8, hjust = 0, color = "red") + # nudge asterisk to the right
    facet_wrap(~param, scales = "free_x") +
    labs(x = "Effect size (95% CI)", y = NULL) +
    scale_color_manual(values = c("Significant" = "darkred", "Not Significant" = "gray40")) +
    theme_minimal() +
    theme( 
        text = element_text(family = "Times", colour = "black"),
        legend.title = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        panel.border = element_rect(color = "black", size = 0.5, fill = "transparent")
    )

# save model averaged effect sizes 
ggsave("Global/Figures/SingleSpeciesModeling/ModelAveragedEffectSizes.png", width = 10, height = 6)


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

# Step 6: Print with kable
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
            "Community", "percentNatural", "RainfallScaled",
            "DistToWater", "TemperatureScaled", "DistToComm"
        )
        siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
        siteCovariate$Rainfall <- siteCovariate$Rainfall * 1000
        allCommunities <- unique(siteCovariate$Community)
        N <- 50
        dfTemplate <- data.frame(
            Community = rep(allCommunities, each = N), 
            RainfallScaled = mean(siteCovariate$RainfallScaled),
            percentNatural = mean(siteCovariate$percentNatural),
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
                          summarize(perc = mean(percentNatural)) %>% 
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
peccary <- ~ atop(paste("Collared peccary"), paste("(", italic("Pecari tajacu"), ")"))
brocket <- ~ atop(paste("Brocket"), paste("(", italic("Mazama sp."), ")"))
paca <- ~ atop(paste("Lowland paca"), paste("(", italic("Cuniculus paca"), ")"))
trumpeter <- ~ atop(paste("Grey-winged trumpeter"), paste("(", italic("Psophia crepitans"), ")"))
fourEyed <- ~ atop(paste("Brown four-eyed opossum"), paste("(", italic("Metachirus nudicaudatus"), ")"))
agouti <- ~ atop(paste("Black agouti"), paste("(", italic("Dasyprocta fuliginosa"), ")"))
armadillo <- ~ atop(paste("Nine-banded armadillo"), paste("(", italic("Dasypus novemcinctus"), ")"))
tinamou <- ~ atop(paste("Great tinamou"), paste("(", italic("Tinamus major"), ")"))
opossum <- ~ atop(paste("Common opossum"), paste("(", italic("Didelphis marsupialis"), ")"))
ocelot <- ~ atop(paste("Ocelot"), paste("(", italic("Leopardus pardalis"), ")"))

# rphylopic per species
peccPic <- get_phylopic(uuid = get_uuid(name = "Pecari tajacu", n = 1))
brockPic <- get_phylopic(uuid = get_uuid(name = "Mazama americana", n = 1))
pacaPic <- get_phylopic(uuid = get_uuid(name = "Cuniculus paca", n = 1))
trumpPic <- get_phylopic(uuid = get_uuid(name = "Psophia crepitans", n = 1))
#fourEyedPic <- get_phylopic(uuid = get_uuid(name = "Metachirus nudicaudatus", n = 1))
agoutiPic <- get_phylopic(uuid = get_uuid(name = "Dasyprocta", n = 1))
armadilloPic <- get_phylopic(uuid = get_uuid(name = "Dasypus novemcinctus", n = 1))
tinamouPic <- get_phylopic(uuid = get_uuid(name = "Tinamus major", n = 1))
#opossumPic <- get_phylopic(uuid = get_uuid(name = "Didelphis", n = 1))
ocelotPic <- get_phylopic(uuid = get_uuid(name = "Leopardus pardalis", n = 1))

# plot null occupancy
# dodge <- position_dodge(width = 0.3)
p <- ggplot(estimates, aes(x = Species,
                              y = nullOccupancy)) +
                              #color = Community)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = nullOccupancy - nullOccupancySE, 
                    ymax = nullOccupancy + nullOccupancySE), 
                 width = 0.15, linewidth = .5) +
  #scale_color_manual(values = c("darkorange", "royalblue", "green3", "yellow3")) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter, agouti, armadillo, tinamou, ocelot)) +
  labs(x = "Species", y = "Null occupancy probability (SE)") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.60),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none",
        panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)) + 
  add_phylopic(peccPic, alpha = 0.2, x = 1.0, y = 0.05, ysize = 0.1) +
  add_phylopic(brockPic, alpha = 0.2, x = 2.0, y = 0.05, ysize = 0.125) +
  add_phylopic(pacaPic, alpha = 0.2, x = 3.0, y = 0.05, ysize = 0.1) +
  add_phylopic(trumpPic, alpha = 0.2, x = 4.0, y = 0.05, ysize = 0.13) +
  #add_phylopic(fourEyedPic, alpha = 0.2, x = 5.0, y = 0.05, ysize = 0.1) +
    add_phylopic(agoutiPic, alpha = 0.2, x = 5.0, y = 0.05, ysize = 0.1) +
    add_phylopic(armadilloPic, alpha = 0.2, x = 6.0, y = 0.05, ysize = 0.1) +
    add_phylopic(tinamouPic, alpha = 0.2, x = 7.0, y = 0.05, ysize = 0.125) +
    #add_phylopic(opossumPic, alpha = 0.2, x = 9.0, y = 0.05, ysize = 0.1) +
    add_phylopic(ocelotPic, alpha = 0.2, x = 8.0, y = 0.05, ysize = 0.1)

# plot with the animal silhouettes
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
# You can also specify other variables to include in the model, like this:
# model <- lm(occupancy ~ treatment + other_variable, data = your_data)
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
    add_phylopic(peccPic, alpha = 0.2, x = 1.0, y = 0.05, ysize = 0.1) +
    add_phylopic(brockPic, alpha = 0.2, x = 2.0, y = 0.05, ysize = 0.125) +
    add_phylopic(pacaPic, alpha = 0.2, x = 3.0, y = 0.05, ysize = 0.1) +
    add_phylopic(trumpPic, alpha = 0.2, x = 4.0, y = 0.05, ysize = 0.13) +
    #add_phylopic(fourEyedPic, alpha = 0.2, x = 5.0, y = 0.05, ysize = 0.1) +
    add_phylopic(agoutiPic, alpha = 0.2, x = 5.0, y = 0.05, ysize = 0.1) +
    add_phylopic(armadilloPic, alpha = 0.2, x = 6.0, y = 0.05, ysize = 0.1) +
    add_phylopic(tinamouPic, alpha = 0.2, x = 7.0, y = 0.05, ysize = 0.125) +
    #add_phylopic(opossumPic, alpha = 0.2, x = 9.0, y = 0.05, ysize = 0.1) +
    add_phylopic(ocelotPic, alpha = 0.2, x = 8.0, y = 0.05, ysize = 0.1)

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
  } else if (covariatesMinusCommunity[i] == "percentNatural"){
    xlabel <- "Percent natural area within 25 km"
  } else if (covariatesMinusCommunity[i] == "DistToComm") {
    xlabel <- "Distance to a community (km)"
  } else if (covariatesMinusCommunity[i] == "DaysEffortScaled") {
    xlabel <- "Camera trap effort (days, scaled)"
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





################################################################################
########################### PLOT MAP OF STATIONS ###############################
################################################################################

# load data
data("world")
stations <- read.csv("Global/Data/AllStationsFormatted.csv")
stations <- stations %>%
  select(c(Station, gps_x, gps_y, CommunityName)) %>%
  distinct()
stations$CommunityName <- gsub("Zabalo", "Zábalo", x = stations$CommunityName)

# order the communities by percent of natural cover
stations <- stations %>% dplyr::filter(CommunityName %in% orderedCommunities)
stations$Community <- factor(stations$CommunityName, levels=orderedCommunities)
stations$CommunityName <- NULL

# only highlight Ecuador
SA <- c("ecuador", "bolivia", "brazil", "chile", "colombia", "argentina", "guyana", "paraguay", "peru", "suriname", "uruguay", "venezuela")
mapColors <- rep("white", length(SA))
mapColors[2] <- "lightyellow"
colors <- c(
    "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
)

# option A with data("world")
worldEdited <- world
worldEdited$Ecu <- ifelse(worldEdited$name_long == "Ecuador", "A", "B")

# option B with 'maps' package
mappy <- maps::map("worldHires", SA)
mappy_sf <- st_as_sf(mappy, crs = st_crs(worldEdited), fill = FALSE)
mappy_sf$Ecu <- ifelse(mappy_sf$ID == "Ecuador", "A", "B")

# option A: less detail in Ecuador border
ecuadorMap <- worldEdited %>%
  dplyr::filter(continent == "South America") %>%
  ggplot() +
  geom_sf(aes(fill=Ecu)) +
  geom_point(data = stations, aes(gps_x, gps_y, fill = Community), pch = 21) +
  #scale_fill_manual(values = c("lightyellow", "white")) +
  scale_fill_manual(name = "Community",
                    values = c(colors, "A" = "lightyellow", "B" = "white"), 
                    breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")) +
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(-150, -37)) + 
  #guides(fill = "none") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"))

# option B: more detail in Ecuador border but can't get Ecuador to highlight yellow?
ecuadorMap <- mappy_sf %>%
  ggplot() +
  geom_sf(aes(fill=Ecu), lwd = 0.5) +
  geom_point(data = stations, aes(gps_x, gps_y, fill = Community), pch = 21) +
  #scale_fill_manual(values = c("lightyellow", "white")) +
  scale_fill_manual(name = "Community",
                    values = c(colors, "A" = "lightyellow", "B" = "white"), 
                    breaks = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 6.5, fill = colors) )) +
  coord_sf(default_crs = sf::st_crs(4326), xlim = c(-150, -37)) + 
  #guides(fill = "none") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"))

# set bounding box to magnify
coords <- as.matrix(stations[,c("gps_x","gps_y")])
e <- as.vector(ext(coords)) 
e["xmin"] <- e["xmin"] - 0.1
e["ymin"] <- e["ymin"] - 0.5
e["xmax"] <- e["xmax"] + 0.1
e["ymax"] <- e["ymax"] + 0.5

# plot map inlay
together <- ecuadorMap + geom_magnify(from = e, 
                          to = c(xmin = -150, xmax = -85, ymin = -45, ymax = 3))
together

if (savePlots == "YES") {
  ggsave(plot = together, filename = paste0(communities, "/Figures/mapInlayWithSites.png"), 
         width = 7, height = 5)
}




################################################################################
######################### PLOT PERCENT NATURAL AREA ############################
################################################################################

# plot percent natural area within 25 km with SD
siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
natStats <- siteCovariate %>% 
  group_by(Community) %>%
  summarize(avgNat = mean(percentNatural), sdNat = sd(percentNatural))
natStats$Community <- gsub("Zabalo", "Zábalo", x = natStats$Community)

# order the communities by percent of natural cover
natStats <- natStats %>% dplyr::filter(Community %in% orderedCommunities)
natStats$Community <- factor(natStats$Community, levels=orderedCommunities)

# plot it
ggplot(natStats, aes(x = Community, y = avgNat, fill = Community)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = avgNat - sdNat, ymax = avgNat + sdNat), width = 0.2) +
  ylab("Percent natural area within 25 km") +
  scale_fill_manual(values = colors) +
  ylim(c(0,1)) +
  theme_bw()+
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"))
# save it
if (savePlots == "YES") {
  ggsave(filename = paste0(communities, "/Figures/percentNatArea25km.png"), 
         width = 7, height = 5)
}

# zoomed
ggplot(natStats, aes(x = Community, y = avgNat, fill = Community)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = avgNat - sdNat, ymax = avgNat + sdNat), width = 0.2) +
  ylab("Percent natural area within 25 km") +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(0.825,1)) +
  theme_bw() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"))
# save it
if (savePlots == "YES") {
  ggsave(filename = paste0(communities, "/Figures/percentNatArea25kmZoomed.png"), 
         width = 7, height = 5)
}






################################################################################
####################### SITE ABUNDANCE AND DIVERSITY ###########################
########################### AND CAMERA TRAP INFO ###############################
################################################################################

# load data
Data <- read.csv("Global/Data/AllIndependentRecordsFormatted.csv") 
Traps <- read.csv("Global/Data/AllStationsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
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
         PercentNaturalArea = mean(siteCovariate$percentNatural[siteCovariate$Community == "Sinangoe"]),
         OperatingDays = round(as.numeric(max(noUnknownsSGE$DateTimeOriginal)-
                                            min(noUnknownsSGE$DateTimeOriginal))))
wholeDiversityZAB <- noUnknownsZAB %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Zábalo", 
         PercentNaturalArea = round(mean(siteCovariate$percentNatural[siteCovariate$Community == "Zabalo"]), 3), 
         OperatingDays = round(as.numeric(max(noUnknownsZAB$DateTimeOriginal)-
                                            min(noUnknownsZAB$DateTimeOriginal))))
wholeDiversitySNA <- noUnknownsSNA %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Siona", 
         PercentNaturalArea = round(mean(siteCovariate$percentNatural[siteCovariate$Community == "Siona"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsSNA$DateTimeOriginal)-
                                            min(noUnknownsSNA$DateTimeOriginal))))
wholeDiversitySPA <- noUnknownsSPA %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "San Pablo", 
         PercentNaturalArea = round(mean(siteCovariate$percentNatural[siteCovariate$Community == "San Pablo"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsSPA$DateTimeOriginal)-
                                            min(noUnknownsSPA$DateTimeOriginal))))

wholeDiversityREM <- noUnknownsREM %>%
  group_by(Species) %>%
  summarise(abundance = n()) %>%
  mutate(Community = "Remolino", 
         PercentNaturalArea = round(mean(siteCovariate$percentNatural[siteCovariate$Community == "Remolino"]), 3),  
         OperatingDays = round(as.numeric(max(noUnknownsREM$DateTimeOriginal)-
                                            min(noUnknownsREM$DateTimeOriginal))))

# abundance and diversity for all communities
communityAbundance <- rbind(wholeDiversityZAB, wholeDiversitySPA, wholeDiversityREM, 
                            wholeDiversitySGE, wholeDiversitySNA)

communityDiversity <- communityAbundance %>%
  group_by(Community, PercentNaturalArea) %>%
  summarise(nIndiv=sum(abundance),
            nSpecies = length(unique(Species)),
            OperatingDays = mean(OperatingDays),
            shannonIndex = round(-sum((abundance/sum(abundance))*log(abundance/sum(abundance))), 3),
            simpsonIndex = round(1-sum((abundance/sum(abundance))^2), 3)) 
communityDiversity$PercentNaturalArea <- round(communityDiversity$PercentNaturalArea, 3)
communityDiversity <- arrange(communityDiversity, desc(PercentNaturalArea))
communityDiversity$Community <- factor(communityDiversity$Community, 
                                       levels = communityDiversity$Community)
communityDiversity
write.csv(communityDiversity, "Global/Data/CommunityDiversityAbundance.csv")



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
  kbl(communityDiversity, col.names = c("Community", "Proportion of Natural Area", 
                                        "Number of Detections", "Number of Species",
                                        "Number of Sampling Days",
                                        "Shannon Diversity Index", "Simpson Diversity Index")) %>%
    kable_classic(full_width = T, html_font = "TimesNewRoman") %>%
    kableExtra::save_kable(file = "Global/Figures/communityDiversityAbundance.png", zoom = 1.5)
  
  # table with just diversity information
  kbl(communityDiversity[,c("Community", "PercentNaturalArea", "OperatingDays", 
                            "shannonIndex", "simpsonIndex")], 
      col.names = c("Community", "Proportion of Natural Area", "Number of Sampling Days", 
                    "Shannon Diversity Index", "Simpson Diversity Index")) %>%
    kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
    save_kable(file = "Global/Figures/communityDiversitySummary.png", zoom = 2)
  
  # table with camera trap information
  kbl(cameraInfo, col.names = c("Community", "Number of Sampling Days", 
                                "Sampling Start Date", "Sampling End Date",
                                "Number of Sites")) %>%
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
      mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCut, col.names = c("Species", "Days per Detection Occasion",
                                 "Average Occupancy Estimate", "Occupancy SE", 
                                 "Average Detection Estimate", "Detection SE")) %>%
      kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
      kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/occupancyDetectionEstimates.png", zoom = 10)
    
    estimatesCutNull <- estimatesWithDays |>
        dplyr::select(CommonNames, DaysGrouped, nullOccupancy, nullOccupancySE, nullDetection, nullDetectionSE) |>
        arrange(desc(nullOccupancy)) |>
        mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCutNull, col.names = c(
        "Species", "Days per Detection Occasion",
        "Null Occupancy Estimate", "Null Occupancy SE",
        "Null Detection Estimate", "Null Detection SE"
    )) %>%
        kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
        kableExtra::save_kable(file = "Global/Figures/SingleSpeciesModeling/nullOccupancyDetectionEstimates.png", zoom = 10)

}


# Make a table with the covariates included in the top models for each species
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



# TIME!
#toc() 


