
################################################################################
############# OCCUPANCY ESTIMATES ACROSS COMMUNITIES AND SPECIES ###############
################################################################################

# mission: 
# requirements: 

setwd("~/Documents/amazon")


# load packages
require(unmarked)
require(dplyr)
require(ggplot2)
require(reshape2)


# input
communities <- c("Sinangoe", "Siona", "Siekopai", "Zabalo")
communitiesAbrv <- c("SGE", "SNA", "SKP", "ZAB")
  # Sinangoe = SGE
  # Siona = SNA
  # Siekopai = SKP
  # Zabalo = ZAB
speciesNames <- c("Cuniculus paca", "Mazama americana", "Dicotyles tajacu", "Psophia crepitans")
commonNames <- c("Collared peccary", "Red brocket", "Lowland paca", "Grey-winged trumpeter") # listTitles
  # paca = Cuniculus paca
  # brocket = Mazama americana
  # collared peccary = Dicotyles tajacu 
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
    siteCovariate <- data.frame(DistToComm = scale(stations$Community/1000)) # site covariates (scaled)
      } else if (communities[i] == "Zabalo") {
        load('Zabalo/Data/R Objects/siteCovs2018.RData') # loads 'siteCovariate'
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

# REMAINING:
  # detection: c("1")
  # occupancy: c("1)

# output lists:
masterBestofTheBest <- list() # occu output for the best model for each species
masterTopModels <- list() # df of the models within 2 AIC of the top model
masterBestModsFitLists <- list() # fit list of the best models for model averaging

for (i in 1:length(communities)) {
  
  # set the variables to be considered for America's Next Top Model
  if (communities[i] == "Zabalo") {
    match_detection <- c("Effort", "Habitat")
    match_occupancy <- c("Habitat", "Trail.Distance", "HuntingIntensity")
  } else if (communities[i] == "Sinangoe"){
    match_detection <- c("1")
    match_occupancy <- c("DistToComm")
  } else {
    match_detection <- c("1")
    match_occupancy <- c("1")
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
  }
  names(bestModsFitLists) <- speciesNames
  
  masterBestModsFitLists[[i]] <- bestModsFitLists
  
  # just the #1 model for each species
  bestOfTheBest <- list()
  for (j in 1:length(speciesNames)) {
    bestOfTheBest[[j]] <- occu(formula(topModels[[j]]$ModelName[1]), ufoMasterList[[i]][[j]])
  }
  names(bestOfTheBest) <- speciesNames
  
  masterBestofTheBest[[i]] <- bestOfTheBest
  
}


  
################################################################################
########################## ESTIMATE CALCULATING ################################
################################################################################

# output lists:
masterEstimatedParameters <- list() # df of occ/det and their SE/range
masterUnmarkedPredOcc <- list() # unmarked::predict() output
masterUnmarkedPredDet <- list() # unmarked::predict() output

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
    
  } else {
    
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

    
  }
  names(estimatedParameters) <- speciesNames
  

  masterEstimatedParameters[[i]] <- estimatedParameters
  masterUnmarkedPredOcc[[i]] <- unmarkedPredOcc
  masterUnmarkedPredDet[[i]] <- unmarkedPredDet
  
} 
  
names(masterEstimatedParameters) <- communities
names(masterUnmarkedPredOcc) <- communities
names(masterUnmarkedPredDet) <- communities




