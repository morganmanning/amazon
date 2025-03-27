
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
require(rphylopic)


# input
communities <- c("Sinangoe", "Zabalo", "Remolino", "San Pablo", "Siona")
communitiesAbrv <- c("SGE", "ZAB", "REM", "SPA", "SNA")
  # Sinangoe = SGE
  # Siona = SNA
  # Siekopai = SKP
  # Zabalo = ZAB
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

    for (j in 1:length(speciesNames)) {
        csv <- paste0(communitiesAbrv[i], gsub(" ", "", speciesNames[j]), ".csv")
        species <- read.csv(paste0(communities[i], "/Data/", csv))
        species <- species[, -1] # remove the column of ID names

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
    if (communities[i] == "Sinangoe") {
        siteCovariate <- data.frame(DistToComm = c(scale(stations$Community / 1000))) # scale
    } else if (communities[i] == "Zabalo") {
        load("Zabalo/Data/R Objects/siteCovs2018.RData") # loads 'siteCovariate'
    } else if (communities[i] == "Global") {
        siteCovariate <- read.csv("Global/Data/AllCommunityCovariates.csv")
        siteCovariate$Rainfall <- siteCovariate$Rainfall * 1000 # convert to grams/m^2/s
        siteCovariate$Community <- factor(siteCovariate$Community,
            levels = c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
        )
    } else {
        siteCovariate <- NULL # no covariates for remaining communities
    }

    # turn all species detection histories for community i into matrices
    for (j in 1:length(detHistory)) {
        detHistory[[j]] <- detHistory[[j]][order(as.numeric(row.names(detHistory[[j]]))), ] # order matters
        detHistory[[j]] <- as.matrix(detHistory[[j]])
    }

    # clump all matrices according to their best clumping factor
    source("./Zabalo/Scripts/optimalClumping.R")

    clumpedMatrixList <- list()
    for (j in 1:length(detHistory)) { # for each species
        y <- detHistory[[j]] # detection history for each species
        clumpEvery <- as.numeric(best_clumping_factor(y)[1])
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
        clumpedMatrixList[[j]] <- clumpedMatrix
    }
    names(clumpedMatrixList) <- commonNames # now we have a matrix of clumped detection histories

    # make unmarked df for each species with clumping
    ufoList <- list()
    for (j in 1:length(clumpedMatrixList)) {
        clumpedMatrix <- clumpedMatrixList[[j]]
        ufo <- unmarkedFrameOccu(clumpedMatrix,
            siteCovs = siteCovariate,
            obsCovs = NULL
        )
        ufoList[[j]] <- ufo
    }
    names(ufoList) <- commonNames

    # make unmarked df for each species without clumping
    unclumpedUFOList <- list()
    for (j in 1:length(detHistory)) {
        unclumpedMatrix <- detHistory[[j]]
        ufo <- unmarkedFrameOccu(unclumpedMatrix,
            siteCovs = siteCovariate,
            obsCovs = NULL
        )
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
            scale_fill_manual(values = c("gray75", "red3"), na.value = "white", name = "") +
            scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
            scale_alpha_identity(guide = "none") +
            coord_equal(expand = 0) +
            xlab(paste(
                "Time (1 unit = ~",
                as.numeric(best_clumping_factor(detHistory[[j]])) * 2,
                # number of columns/time steps divided by the clumping factor, times 2
                "days)"
            )) +
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

        # plot
        unclumpPlotList[[j]] <- ggplot(meltedComb, aes(Var2, Var1, fill = value)) +
            geom_tile(colour = "gray50") +
            scale_fill_manual(values = c("gray75", "red3"), na.value = "white", name = "") +
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
  match_detection <- c("1")
  match_occupancy <- c("1")
  
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
    # for no covariates
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
allCommunities <- c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
communitiesAccent <- gsub(pattern = "Zabalo", replacement = "Zábalo", communities)
speciesNames

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
    row <- (estimates$Community == communitiesAccent[i]) & (estimates$Species == speciesNames[j])
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
estimates$Species <- factor(estimates$Species, levels = speciesNames) # so plotting doesn't alphabetize species







################################################################################
################################################################################
################## THINGS ARE NOT AUTOMATED FROM HERE!!! #######################
######################## REQUIRES MANUAL INPUT!!! ##############################
################################################################################
################################################################################

########## MANUALLY MAKE X-AXIS LABELS FOR EACH SPECIES: ENSURE IN CORRECT ORDER 
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
agoutiPic <- get_phylopic(uuid = get_uuid(name = "Dasyprocta", n = 1))
armadilloPic <- get_phylopic(uuid = get_uuid(name = "Dasypus novemcinctus", n = 1))
tinamouPic <- get_phylopic(uuid = get_uuid(name = "Tinamus major", n = 1))
opossumPic <- get_phylopic(uuid = get_uuid(name = "Didelphis", n = 1))
ocelotPic <- get_phylopic(uuid = get_uuid(name = "Leopardus pardalis", n = 1))

#### in the correct order???
commonNames

########## MANUALLY INPUT X-AXIS LABELS IN THE CORRECT ORDER
# plot it
estimates$Community <- factor(estimates$Community,
                                        levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona"))
colors <- c("Zábalo" = "darkgreen", "Remolino" = "forestgreen", 
            "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3")
dodge <- position_dodge(width = 0.3)
plot <- ggplot(estimates, aes(x = Species,
                           y = avgOccupancy,
                           color = Community)) +
  geom_point(aes(color = Community), position = dodge, size = 2.5) +
  geom_errorbar(aes(ymin = avgOccupancy - avgOccupancySE, 
                    ymax = avgOccupancy + avgOccupancySE, 
                    color = Community), 
                position = dodge, width = 0.2, linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(peccary, brocket, paca, trumpeter, fourEyed, agouti, armadillo, tinamou, opossum, ocelot)) +
  labs(x = "Species", y = "Naive occupancy probability") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position="top")
plot

commonNames
# plot the animal silhouettes
plot +
    add_phylopic(peccPic, alpha = 0.2, x = 1.0, y = 0.05, ysize = 0.1) +
    add_phylopic(brockPic, alpha = 0.2, x = 2.0, y = 0.05, ysize = 0.125) +
    add_phylopic(pacaPic, alpha = 0.2, x = 3.0, y = 0.05, ysize = 0.1) +
    add_phylopic(trumpPic, alpha = 0.2, x = 4.0, y = 0.05, ysize = 0.13) +
    add_phylopic(fourEyedPic, alpha = 0.2, x = 5.0, y = 0.05, ysize = 0.1) +
    add_phylopic(agoutiPic, alpha = 0.2, x = 6.0, y = 0.05, ysize = 0.1) +
    add_phylopic(armadilloPic, alpha = 0.2, x = 7.0, y = 0.05, ysize = 0.1) +
    add_phylopic(tinamouPic, alpha = 0.2, x = 8.0, y = 0.05, ysize = 0.125) +
    add_phylopic(opossumPic, alpha = 0.2, x = 9.0, y = 0.05, ysize = 0.1) +
    add_phylopic(ocelotPic, alpha = 0.2, x = 10.0, y = 0.05, ysize = 0.1)

# save it
ggsave(filename = "Global/Figures/AllCommunitiesOccupancyEstimatesNullModels.png", width = 8, height = 4)



