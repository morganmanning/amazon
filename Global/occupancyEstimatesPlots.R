
################################################################################
############# OCCUPANCY ESTIMATES ACROSS COMMUNITIES AND SPECIES ###############
################################################################################
# mission: plot estimates for each species in each community

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

tic() # time it

# input
communities <- "Global"
communitiesAbrv <- "All"
#speciesNames <- c("Pecari tajacu", "Mazama americana", "Cuniculus paca", "Psophia crepitans")
#commonNames <- c("Collared peccary", "Red brocket", "Lowland paca", "Grey-winged trumpeter") # listTitles
speciesNames <- c("Pecari tajacu", "Mazama sp.", "Cuniculus paca", "Psophia crepitans", "Metachirus nudicaudatus", "Dasyprocta fuliginosa", "Dasypus novemcinctus", "Tinamus major", "Didelphis marsupialis", "Leopardus pardalis")
commonNames <- c("Collared peccary", "Brockets", "Lowland paca", "Grey-winged trumpeter", "Brown four-eyed opossum", "Black agouti", "Nine-banded armadillo", "Great tinamou", "Common opossum", "Ocelot") 
  # paca = Cuniculus paca
  # brocket = Mazama americana
  # collared peccary = Pecari tajacu 
  # trumpeter = Psophia crepitans
  # brown four-eyed possum = Metachirus nudicaudatus (#1 species in SGE)
  # black agouti = Dasyprocta fuliginosa (#2 species in SGE)



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
  for (j in 1:length(detHistory)) { # for each species
    y <- detHistory[[j]] # detection history for each species
    clumpEvery <- as.numeric(best_clumping_factor(y)[1])
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



################################################################################
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
bestModsOutputs <- list()
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
    match_detection <- c("Community")
    match_occupancy <- c("Community", "RainfallScaled", "percentNatural",
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
      df[m,1]<- as.character(c(allModels[[j]][[m]]@formula))
      df[m,2]<- allModels[[j]][[m]]@AIC
    }
    df <- df[order(df$AIC),]
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
      test <- occu(~Community ~RainfallScaled + percentNatural + DistToComm + TemperatureScaled + Community, 
                    ufoMasterList[[1]][[j]]) # i = 1 when all communities are together
      
      globalModels[[j]] <- occu(~Community ~RainfallScaled + percentNatural + DistToComm +
                    TemperatureScaled + Community, 
                    ufoMasterList[[1]][[j]], 
                                   control = 10000, 
                                   starts = c(rep(0, length(test@opt$par)))) 
      #globalModelsFitList[[j]] <- fitList(global = globalModels[[j]])
}
names(globalModels) <- speciesNames
globalModelsNested <- list(Global = globalModels) # need to nest so it matches the dimensions of masterBestModsFitLists moving forwards

# to proceed using the global models rather than the best models, uncomment the following
#masterBestModsFitLists <- list()
#masterBestModsFitLists <- globalModelsNested
#names(masterBestModsFitLists) <- communities



  
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
  if (communities[i] == "Sinangoe"){
    siteCovariate <- data.frame(DistToComm = scale(stations$Community/1000)) # site covariates (scaled)
    df <- data.frame(
      DistToComm = siteCovariate$DistToComm
    ) 
    
    for (j in 1:length(speciesNames)) { # for each species
      # state = occupancy
      unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                type = 'state', new = df, appendData = TRUE)
      # det = detection
      bestModel <- occu(formula(masterTopModels[[i]][[j]]$ModelName[1]), 
                        ufoMasterList[[i]][[j]])
      nullDetPred <- backTransform(bestModel, "det")
      
      # average
      estimatedParameters[[j]] <- data.frame(avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                                             avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                                             occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                                             avgDetection = nullDetPred@estimate,
                                             avgDetectionSE = SE(nullDetPred),
                                             detectionRange = NA)  
    }
    
  } else if (communities[i] == "Zabalo") {
    load('Zabalo/Data/R Objects/siteCovs2018.RData') # loads 'siteCovariate'
    df <- data.frame(
      Habitat = siteCovariate$Habitat,
      HuntingIntensity = siteCovariate$HuntingIntensity,
      Trail.Distance = siteCovariate$Trail.Distance,
      Effort = siteCovariate$Effort
    ) 
    
    for (j in 1:length(speciesNames)) { # for each species
      # state = occupancy
      unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                type = 'state', new = df, appendData = TRUE)
      # det = detection
      unmarkedPredDet[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                type = 'det', new = df, appendData = TRUE)
      
      # average
      estimatedParameters[[j]] <- data.frame(avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                                             avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                                             occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                                             avgDetection = mean(unmarkedPredDet[[j]]$Predicted),
                                             avgDetectionSE = mean(unmarkedPredDet[[j]]$SE),
                                             detectionRange = range(unmarkedPredDet[[j]]$Predicted))  
    }
    
  } else if (communities[i] == "Siona" | communities[i] == "Siekopai") {
    
    # for Siona and Siekopai who don't have covariates
    siteCovariate <- NULL # no covariates for remaining communities 
    df <- NULL # not sure what to do with this since we don't have covariates for other communities
    
    # calculate occupancy and detection estimates
    for (j in 1:length(speciesNames)){
      bestModel <- occu(~1 ~1, ufoMasterList[[i]][[j]])
      nullOccPred <- backTransform(bestModel, "state")
      nullDetPred <- backTransform(bestModel, "det")
      
      # average
      estimatedParameters[[j]] <- data.frame(avgOccupancy = nullOccPred@estimate,
                                             avgOccupancySE = SE(nullOccPred),
                                             occupancyRange = NA,
                                             avgDetection = nullDetPred@estimate,
                                             avgDetectionSE = SE(nullDetPred),
                                             detectionRange = NA)  
    }
    
    
  } else if (communities[i] == "Global"){
    covariates <- c("Community", "RainfallScaled", "percentNatural", 
                    "DistToWater", "TemperatureScaled", "DistToComm")
    siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
    siteCovariate$Rainfall <- siteCovariate$Rainfall * 1000
    allCommunities <- unique(siteCovariate$Community)
    N <- 50
    dfTemplate <- data.frame(
      Community = rep("Sinangoe", N), # picked because it's kind of a middle community
      RainfallScaled = mean(siteCovariate$RainfallScaled),
      percentNatural = mean(siteCovariate$percentNatural),
      DistToWater = mean(siteCovariate$DistToWater),
      TemperatureScaled = mean(siteCovariate$TemperatureScaled),
      DistToComm = mean(siteCovariate$DistToComm)
    ) 
    
    for (j in 1:length(speciesNames)) { # for each species
      dfEdited <- dfTemplate
      ########## average overall
      # state = occupancy
      unmarkedPredOcc[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                type = 'state', new = dfTemplate, appendData = TRUE)
      # det = detection
      unmarkedPredDet[[j]] <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                type = 'det', new = dfTemplate, appendData = TRUE)
      
      # average
      estimatedParameters[[j]] <- data.frame(avgOccupancy = mean(unmarkedPredOcc[[j]]$Predicted),
                                             avgOccupancySE = mean(unmarkedPredOcc[[j]]$SE),
                                             occupancyRange = range(unmarkedPredOcc[[j]]$Predicted),
                                             avgDetection = mean(unmarkedPredDet[[j]]$Predicted),
                                             avgDetectionSE = mean(unmarkedPredDet[[j]]$SE),
                                             detectionRange = range(unmarkedPredDet[[j]]$Predicted))
      for (m in 1:length(allCommunities)) {
        dfEdited <- dfTemplate
        dfEdited[,"Community"] <- allCommunities[m]
        
        ######## prediction per covariate
        for (k in 1: length(covariates)) {
          dfEdited <- dfTemplate
          dfEdited[,"Community"] <- allCommunities[m]
          covariateInQuestion <- covariates[k]
          
          # prediction DF
          if (covariateInQuestion == "Community"){
            dfEdited[,covariateInQuestion] <- rep(unique(siteCovariate$Community), 
                                                  each = N/length(unique(siteCovariate$Community)))
          } else {
            dfEdited[,"Community"] <- allCommunities[m]
            dfEdited[,covariateInQuestion] <- seq(min(siteCovariate[,covariateInQuestion]), 
                                                  max(siteCovariate[,covariateInQuestion]), 
                                                  length.out = N)
          }
          
          # state = occupancy
          df <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                    type = 'state', new = dfEdited, appendData = TRUE)
          df$PredictedCovariate <- covariateInQuestion
          df$CommunityHeldConstant <- allCommunities[m]
          df$Species <- speciesNames[j]
          perCovOccupancy[[k]] <- df
          
          # det = detection
          df <- unmarked::predict(masterBestModsFitLists[[i]][[j]], 
                                                    type = 'det', new = dfEdited, appendData = TRUE)
          df$PredictedCovariate <- covariateInQuestion
          df$CommunityHeldConstant <- allCommunities[m]
          df$Species <- speciesNames[j]
          perCovDetection[[k]] <- df
          
        }
        names(perCovOccupancy) <- covariates
        names(perCovDetection) <- covariates
        
        perCovPerCommOccupancy[[m]] <- perCovOccupancy
        perCovPerCommDetection[[m]] <- perCovDetection
      }
      names(perCovPerCommOccupancy) <- allCommunities
      names(perCovPerCommDetection) <- allCommunities
      
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
    #estimates[row, "avgOccupancySEManual"] <- sd(masterUnmarkedPredOcc[[i]][[j]]$Predicted)/sqrt(nrow(masterUnmarkedPredOcc[[i]][[j]]))
    #estimates[row, "avgOccupancySD"] <- sd(masterUnmarkedPredOcc[[i]][[j]]$Predicted)
    estimates[row, "avgDetection"] <- masterEstimatedParameters[[i]][[j]]$avgDetection[1]
    estimates[row, "avgDetectionSE"] <- masterEstimatedParameters[[i]][[j]]$avgDetectionSE[1]
    #estimates[row, "avgDetectionSEManual"] <- sd(masterUnmarkedPredDet[[i]][[j]]$Predicted)/sqrt(nrow(masterUnmarkedPredDet[[i]][[j]]))
    #estimates[row, "avgDetectionSD"] <- sd(masterUnmarkedPredDet[[i]][[j]]$Predicted)
  }
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




# FROM HERE: PLOTS WILL ONLY WORK IF SPECIES == PECCARY, BROCKET, PACA, TRUMPETER
commonNames





################################################################################
############################## PLOT ESTIMATES ##################################
################################################################################


# make a dataframe that has predictions across each covariate to facilitate plotting
covariatesMinusCommunity <- covariates[covariates != "Community"] # everything but community

# making the data frame to plot from
plottingDF <- as.data.frame(do.call(rbind, do.call(rbind, do.call(rbind, masterGlobalOccupancyEstimates))))
plottingDF <- plottingDF[plottingDF$PredictedCovariate != "Community",]
plottingDF$Community <- gsub("Zabalo", "Zábalo", x = plottingDF$Community)

# order the communities by percent of natural cover
orderedCommunities <- c(siteCovariate %>% 
                          group_by(Community) %>% 
                          summarize(perc = mean(percentNatural)) %>% 
                          arrange(desc(perc)) %>% 
                          select(Community))$Community
orderedCommunities <- gsub("Zabalo", "Zábalo", x = orderedCommunities)
plottingDF <- plottingDF %>% dplyr::filter(Community %in% orderedCommunities)
plottingDF$Community <- factor(plottingDF$Community, levels=orderedCommunities)

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
fourEyedPic <- get_phylopic(uuid = get_uuid(name = "Metachirus nudicaudatus", n = 1))
agoutiPic <- get_phylopic(uuid = get_uuid(name = "Dasyprocta fuliginosa", n = 1))
armadilloPic <- get_phylopic(uuid = get_uuid(name = "Dasypus novemcinctus", n = 1))
tinamouPic <- get_phylopic(uuid = get_uuid(name = "Tinamus major", n = 1))
opossumPic <- get_phylopic(uuid = get_uuid(name = "Didelphis marsupialis", n = 1))
ocelotPic <- get_phylopic(uuid = get_uuid(name = "Leopardus pardalis", n = 1))

# plot it
dodge <- position_dodge(width = 0.3)
plot <- ggplot(estimates, aes(x = Species,
                              y = avgOccupancy,
                              color = Community)) +
  geom_point(aes(color = Community), position = dodge, size = 2.5) +
  geom_errorbar(aes(ymin = avgOccupancy - avgOccupancySE, 
                    ymax = avgOccupancy + avgOccupancySE, 
                    color = Community), 
                position = dodge, width = 0.2, linewidth = 1) +
  #scale_color_manual(values = c("darkorange", "royalblue", "green3", "yellow3")) +
  scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter, fourEyed, agouti, armadillo, tinamou, opossum, ocelot)) +
  labs(x = "Species", y = "Naive occupancy probability estimate") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position="none") + 
  add_phylopic(peccPic, alpha = 0.2, x = 1.0, y = 0.10, ysize = 0.25) +
  add_phylopic(brockPic, alpha = 0.2, x = 2.0, y = 0.19, ysize = 0.45) +
  add_phylopic(pacaPic, alpha = 0.2, x = 3.0, y = 0.10, ysize = 0.25) +
  add_phylopic(trumpPic, alpha = 0.2, x = 4.0, y = 0.19, ysize = 0.45)

# plot with the animal silhouettes
plot 

# save it
if (savePlots == "YES") {
  ggsave(filename = paste0(communities, "/Figures/", 
                           communitiesAbrv, "OccupancyEstimates.png"), 
         width = 8, height = 4)
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
          #legend.position="top"
    )
  #theme(plot.title = element_text(hjust = 0.5))
  
  # save it
  if (savePlots == "YES") {
    ggsave(filename = paste0(communities, "/Figures/", 
                             cov, "Prediction.png"), 
           width = 8, height = 4)
  }
  
}






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
  ggsave(filename = paste0(communities, "/Figures/mapInlayWithSites.png"), 
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
Data$Species <- gsub("Mazama americana", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama nemorivaga", "Mazama sp.", Data$Species)
Data$Species <- gsub("Mazama gouazoubira", "Mazama sp.", Data$Species) # replace all Mazama species with Mazama sp.


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
    estimatesCut <- estimates |> 
      dplyr::select(CommonNames, avgOccupancy, avgOccupancySE, avgDetection, avgDetectionSE) |>
      arrange(desc(avgOccupancy)) |>
      mutate(across(where(is.numeric), round, 3))
    kbl(estimatesCut, col.names = c("Species",
                                 "Average Occupancy Estimate", "Occupancy SE", 
                                 "Average Detection Estimate", "Detection SE")) %>%
      kable_classic(full_width = FALSE, html_font = "TimesNewRoman") %>%
      kableExtra::save_kable(file = "Global/Figures/occupancyDetectionEstimates.png", zoom = 10)

}


# Make a table with the covariates included in the top models for each species
if (communities == "Global" & savePlots == "YES" & all(speciesNames == names(masterTopModels[[1]])) & length(speciesNames) == 4){
    
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
    nModelsPerSpecies <- unique(bestModelsDF$nModels)
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
            start_row = nrow(bestModelsDF) - nModelsPerSpecies[4] + 1,
            end_row = nrow(bestModelsDF)
        ) %>%
        kableExtra::save_kable(file = "Global/Figures/AllSpeciesBestModelsTable.png", zoom = 10)


}



# TIME!
toc() # typically takes 210 seconds (3.5 minutes)


