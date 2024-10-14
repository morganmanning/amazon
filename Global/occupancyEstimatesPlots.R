
################################################################################
############# OCCUPANCY ESTIMATES ACROSS COMMUNITIES AND SPECIES ###############
################################################################################

# mission: plot estimates for each species in each community
# requirements: Siona and Siekopai have no covariates at this point

setwd("~/Documents/amazon")


# load packages
require(unmarked)
require(dplyr)
require(ggplot2)
require(reshape2)
require(rphylopic)
require(ggpubr)


# input
#communities <- c("Sinangoe", "Siona", "Siekopai", "Zabalo")
#communitiesAbrv <- c("SGE", "SNA", "SKP", "ZAB")
  # Sinangoe = SGE
  # Siona = SNA
  # Siekopai = SKP
  # Zabalo = ZAB
communities <- "Global"
communitiesAbrv <- "All"
speciesNames <- c("Pecari tajacu", "Mazama americana", "Cuniculus paca", "Psophia crepitans")
commonNames <- c("Collared peccary", "Red brocket", "Lowland paca", "Grey-winged trumpeter") # listTitles
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
  
  # import site covariates for each community
  if (communities[i] == "Sinangoe"){
    siteCovariate <- data.frame(DistToComm = c(scale(stations$Community/1000))) # scale
      } else if (communities[i] == "Zabalo") {
        load('Zabalo/Data/R Objects/siteCovs2018.RData') # loads 'siteCovariate'
        } else if (communities[i] == "Global") {
         siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
         siteCovariate$Rainfall <- siteCovariate$Rainfall*1000 # convert to grams/m^2/s
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
    covariates <- c("Community", "RainfallScaled", "percentNatural", "DistToWater", "TemperatureScaled", "DistToComm")
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
################## THINGS ARE NOT AUTOMATED FROM HERE!!! #######################
######################## REQUIRES MANUAL INPUT!!! ##############################
################################################################################
################################################################################

########## MANUALLY MAKE X-AXIS LABELS FOR EACH SPECIES: ENSURE IN CORRECT ORDER 
# work around to make x-axis labels with italics and divided into two lines
speciesNames
peccary <- ~ atop(paste("Collared peccary"), paste("(", italic("Pecari tajacu"), ")"))
brocket <- ~ atop(paste("Red brocket"), paste("(", italic("Mazama americana"), ")"))
paca <- ~ atop(paste("Lowland paca"), paste("(", italic("Cuniculus paca"), ")"))
trumpeter <- ~ atop(paste("Grey-winged trumpeter"), paste("(", italic("Psophia crepitans"), ")"))

#### in the correct order???
speciesNames

########## MANUALLY INPUT X-AXIS LABELS IN THE CORRECT ORDER
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
  scale_color_manual(values = c("darkorange", "royalblue", "green3", "yellow3")) +
  scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter)) +
  labs(x = "Species", y = "Probability of Occupancy") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position="top")
#plot

# add animal silhouettes 
peccPic <- get_phylopic(uuid = get_uuid(name = "Dicotyles tajacu", n = 1))
brockPic <- get_phylopic(uuid = get_uuid(name = "Mazama americana", n = 1))
pacaPic <- get_phylopic(uuid = get_uuid(name = "Cuniculus paca", n = 1))
trumpPic <- get_phylopic(uuid = get_uuid(name = "Psophia crepitans", n = 1))

# plot the animal silhouettes
plot + 
  add_phylopic(peccPic, alpha = 0.2, x = 1.0, y = 0.10, ysize = 0.25) +
  add_phylopic(brockPic, alpha = 0.2, x = 2.0, y = 0.19, ysize = 0.45) +
  add_phylopic(pacaPic, alpha = 0.2, x = 3.0, y = 0.10, ysize = 0.25) +
  add_phylopic(trumpPic, alpha = 0.2, x = 4.0, y = 0.19, ysize = 0.45)

# save it
ggsave(filename = "Global/Figures/GlobalOccupancyEstimates.png", width = 8, height = 4)



######################### GLOBAL 
require(ggpubr)
# GOAL:
# - one plot for each covariate with four panels (one for each species) with prediction lines for each community
covariatesMinusCommunity <- c("RainfallScaled", "percentNatural", "DistToWater", "TemperatureScaled", "DistToComm")

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
covPlots <- list()
allCovPlots <- list()
for (i in 1:length(covariatesMinusCommunity)){
  
  # x axis labels for each covariate
  if (covariatesMinusCommunity[i] == "RainfallScaled") {
    xlabel <- expression(paste0("Rainfall (kg*", m^{-2}, "*", s^{-1}, ", scaled)"))
  } else if (covariatesMinusCommunity[i] == "TemperatureScaled"){
    xlabel <- "Temperature (°C, scaled)"
  } else if (covariatesMinusCommunity[i] == "DistToWater") {
    xlabel <- "Distance to water source (m)"
  } else if (covariatesMinusCommunity[i] == "percentNatural"){
    xlabel <- "Percent natural area within 25 km"
  } else if (covariatesMinusCommunity[i] == "DistToComm") {
    xlabel <- "Distance to a community (km)"
  }
  
  for (j in 1:length(speciesNames)) {
    # subset it to just the species j
    perCovSpp <- plottingDF[plottingDF$Species == speciesNames[j] & 
                              plottingDF$PredictedCovariate == covariatesMinusCommunity[i],]
    
    # get phylopic and look at just the covariate in question
    animalPic <- get_phylopic(uuid = get_uuid(name = speciesNames[j], n = 1))
    covariateInQuestion <- perCovSpp[,covariatesMinusCommunity[i]]
    bottomRightX <- max(covariateInQuestion)-(0.02 * max(covariateInQuestion))
    bottomY <- 0.15
    
    # actually plot it
    plotty <- ggplot(perCovSpp, aes(x = covariateInQuestion, 
                          y = Predicted, color = Community)) +
      geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), 
                  alpha = 0.2, color = NA) +
      geom_line(aes(x = covariateInQuestion, y = Predicted)) +
      add_phylopic(animalPic, alpha = 0.2, x = bottomRightX, y = bottomY, ysize = 0.3) +
      ylab(expression(paste("Occupancy probability estimate (", psi, ")"))) +
      xlab(xlabel) +
      ggtitle(commonNames[j]) +
      coord_cartesian(ylim = c(0,1), 
                      xlim = c(min(covariateInQuestion)-(0.001 * min(covariateInQuestion)),
                               max(covariateInQuestion)+(0.005 * max(covariateInQuestion)))) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    covPlots[[j]] <- plotty
    
  }
  
  names(covPlots) <- commonNames
  
  allCovPlots[[i]] <- covPlots
  
}
names(allCovPlots) <- covariatesMinusCommunity
    
    
    
ggarrange(plotlist = allCovPlots[["Rainfall"]], ncol = 2, nrow = 2, common.legend = TRUE)
    

animalPic <- get_phylopic(uuid = get_uuid(name = speciesNames[j], n = 1))
perCovSpp <- plottingDF[plottingDF$Species == speciesNames[j] & 
                          plottingDF$PredictedCovariate == covariatesMinusCommunity[i],]
covariateInQuestion <- perCovSpp[,covariatesMinusCommunity[i]]
ggplot(perCovSpp, aes(x = covariateInQuestion, 
                   y = Predicted, color = Community)) +
  geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), 
              alpha = 0.2, color = NA) +
  geom_line(aes(x = covariateInQuestion, y = Predicted)) +
  add_phylopic(animalPic, alpha = 0.2, x = 0.98, y = 0.15, ysize = 0.3) +
  ylab(expression(paste("Occupancy probability estimate (", psi, ")"))) +
  xlab(xlabel) +
  coord_cartesian(ylim = c(0,1), 
                  xlim = c(min(covariateInQuestion)-(0.001 * min(covariateInQuestion)),
                                               max(covariateInQuestion)+(0.005 * max(covariateInQuestion)))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_bw()



######## MANUALLY
# animal silhouettes 
peccPic <- get_phylopic(uuid = get_uuid(name = "Dicotyles tajacu", n = 1))
brockPic <- get_phylopic(uuid = get_uuid(name = "Mazama americana", n = 1))
pacaPic <- get_phylopic(uuid = get_uuid(name = "Cuniculus paca", n = 1))
trumpPic <- get_phylopic(uuid = get_uuid(name = "Psophia crepitans", n = 1))



#### RAINFALL 
# just the covariate
df <- plottingDF[plottingDF$PredictedCovariate == "percentNatural",]

covPlots <- list()
for (j in 1:length(speciesNames)) {
  # subset it to just the species j
  perCovSpp <- df[df$Species == speciesNames[j],]
  
  # get phylopic and look at just the covariate in question
  animalPic <- get_phylopic(uuid = get_uuid(name = speciesNames[j], n = 1))
  covariateInQuestion <- perCovSpp[,"percentNatural"]
  bottomRightX <- max(covariateInQuestion)-(0.02 * max(covariateInQuestion))
  bottomY <- 0.15
  
  # actually plot it
  plotty <- ggplot(df, aes(x = percentNatural, 
                                  y = Predicted, color = Community)) +
    geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), alpha = 0.2, color = NA) +
    geom_line(aes(x = percentNatural, y = Predicted)) +
    add_phylopic(animalPic, alpha = 0.2, x = bottomRightX, y = bottomY, ysize = 0.3) +
    ylab(expression(paste("Occupancy probability estimate (", psi, ")"))) +
    #facet_grid(~Species) +
    xlab("Covariate") +
    ggtitle(commonNames[j]) +
    coord_cartesian(ylim = c(0,1), 
                    xlim = c(min(covariateInQuestion),
                             max(covariateInQuestion))) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  covPlots[[j]] <- plotty
  
}


ggarrange(plotlist = covPlots, ncol = 2, nrow = 2, common.legend = TRUE)



##### faceted way
covPlots <- list()
for (i in 1:length(covariatesMinusCommunity)) {
  # subset it to just the species j
  df <- plottingDF[plottingDF$PredictedCovariate == covariatesMinusCommunity[i],]
  covariateInQuestion <- df[,covariatesMinusCommunity[i]]
 
  # actually plot it
  plotty <- ggplot(df, aes(x = covariateInQuestion, 
                           y = Predicted, color = Community)) +
    geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), 
                alpha = 0.2, color = NA) +
    geom_line(aes(x = covariateInQuestion, y = Predicted)) +
    ylab(expression(paste("Occupancy probability estimate (", psi, ")"))) +
    facet_wrap(~commonNames, nrow = 2, ncol = 2) +
    xlab("Covariate") +
    #ggtitle(covariatesMinusCommunity[i]) +
    coord_cartesian(ylim = c(0,1), 
                    xlim = c(min(covariateInQuestion),
                             max(covariateInQuestion))) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() 
    #theme(plot.title = element_text(hjust = 0.5))
  covPlots[[i]] <- plotty
  
}


ggarrange(plotlist = covPlots, ncol = 2, nrow = 2, common.legend = TRUE)



######## manually making and saving prediction plots # this is how I did it for CLAG
# covariate wanted
cov <- "DistToComm" # percentNatural, RainfallScaled, TemperatureScaled, DistToWater, DistToComm

# subset
df <- plottingDF[plottingDF$PredictedCovariate == cov,]
covariateInQuestion <- df[,cov]

# actually plot it
ggplot(df, aes(x = covariateInQuestion, 
               y = Predicted, color = Community)) +
  geom_ribbon(aes(ymin = Predicted - SE, ymax = Predicted + SE, fill = Community), 
              alpha = 0.2, color = NA) +
  geom_line(aes(x = covariateInQuestion, y = Predicted)) +
  ylab(expression(paste("Occupancy probability estimate (", psi, "  )"))) +
  facet_wrap(~commonNames, nrow = 2, ncol = 2) +
  xlab("Distance to a community (km)") +
  coord_cartesian(ylim = c(0,1), 
                  xlim = c(min(covariateInQuestion),
                           max(covariateInQuestion))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_bw() 
#theme(plot.title = element_text(hjust = 0.5))

# save it
ggsave(filename = "Global/Figures/distToComm.png", width = 8, height = 4)


# xlab:
# rainfall: expression(paste("Rainfall (g*", m^{-2}, s^{-1}, ", scaled)"))
# percentNatural: Percent natural area within 25 km
# Temperature: Temperature (°C, scaled)
# DistToWater: Distance to water (m)
# DistToComm: Distance to a community (km)


# plot stations on a map
library(sf)
library(spData)
require(ggmagnify)
require(terra)
require(mapdata)

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
mappy <- maps::map("worldHires", SA, fill = )
mappy_sf <- st_as_sf(mappy, crs = st_crs(worldEdited), fill = FALSE)
mappy_sf$Ecu <- ifelse(mappy_sf$ID == "Ecuador", "A", "B")


worldEdited <- world
worldEdited$Ecu <- ifelse(worldEdited$name_long == "Ecuador", "A", "B")

# make plot
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

# make plot
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
ggsave(plot = together, filename = "Global/Figures/mapInlayWithSites.png", width = 7, height = 5)


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
ggsave("Global/Figures/percentNatArea25km.png", width = 7, height = 5)

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
ggsave("Global/Figures/percentNatArea25kmZoomed.png", width = 7, height = 5)



### hunted species in Zabalo
require(kableExtra)
ZABhunting <- read.csv("../../Zabalo/Data/HuntingData2018.csv")

ZABhuntingSum <- ZABhunting %>%
  group_by(Species) %>%
  summarize(count = sum(Count))

together <- data.frame(Species = c("Peccaries", "Paca", "Trumpeter", "Brockets"),
                       NumberHunted = c(sum(ZABhuntingSum$count[ZABhuntingSum$Species == "White-lipped Peccary" | ZABhuntingSum$Species == "Collared Peccary"]),
                                        sum(ZABhuntingSum$count[ZABhuntingSum$Species == "Lowland Paca"]),
                                        sum(ZABhuntingSum$count[ZABhuntingSum$Species == "Grey-winged Trumpeter"]),
                                        sum(ZABhuntingSum$count[ZABhuntingSum$Species == "Red Brocket Deer" | ZABhuntingSum$Species == "Grey Brocket Deer"])))

kbl(together, 
    col.names = c("Species", "Number Hunted")) %>%
  kable_classic(font_size = 22, html_font = "TimesNewRoman") %>%
  save_kable(file = "../Figures/ZabaloHuntedSpecies.png", zoom = 2)



# trial mapping
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(maps)
library(mapdata)

# Load map data
south_america_map <- map_data("world", region = SA)

# Create a data frame with GPS coordinates of random points in Ecuador
stations <- data.frame(
  Community = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona"),
  gps_x = c(-78.5, -78.2, -79.0, -77.5, -78.0),
  gps_y = c(-1.5, -1.0, -2.0, -3.0, -1.8)
)

# Create a map with Ecuador highlighted
worldEdited <- south_america_map
worldEdited$Ecu <- ifelse(worldEdited$region == "Ecuador", "A", "B")

# Create the main map
ecuadorMap <- ggplot() +
  geom_polygon(data = worldEdited, aes(x = long, y = lat, group = group, fill = Ecu)) +
  geom_point(data = stations, aes(x = gps_x, y = gps_y, fill = Community), pch = 21, size = 4) +
  scale_fill_manual(name = "Community",
                    values = c("A" = "lightyellow", "B" = "white"), 
                    labels = c("Ecuador", "Other Countries")) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"))

# Set bounding box for magnification
coords <- as.matrix(stations[, c("gps_x", "gps_y")])
e <- as.vector(ext(coords))
e["xmin"] <- e["xmin"] - 0.1
e["ymin"] <- e["ymin"] - 0.5
e["xmax"] <- e["xmax"] + 0.1
e["ymax"] <- e["ymax"] + 0.5

# Add a rectangle for zoomed-in area
together <- ecuadorMap + 
  geom_rect(aes(xmin = e["xmin"], xmax = e["xmax"], ymin = e["ymin"], ymax = e["ymax"]), 
            fill = NA, color = "black", linetype = "dashed") +
  annotate("text", x = mean(c(e["xmin"], e["xmax"])), y = mean(c(e["ymin"], e["ymax"])), 
           label = "Zoomed In", size = 5)

# Display the combined plot
print(together)
