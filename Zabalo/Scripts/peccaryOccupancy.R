###### PARAMETER ESTIMATION WORK SHOP ######

setwd("~/Documents/amazon/Data")

require(unmarked)
require(dplyr)
require(ggplot2)

############################################################################
############################# LOAD DATA ####################################
############################################################################

# occupancy data
peccary = read.csv("CollaredPeccary.csv")
peccary = data.frame(peccary[,-1], row.names=peccary[,1]) #rownames = stations

deer = read.csv("Deer.csv")
deer = data.frame(deer[,-1], row.names=deer[,1]) #rownames = stations

paca = read.csv("Paca.csv")
paca = data.frame(paca[,-1], row.names=paca[,1]) #rownames = stations

# occupancy covariates
load('R Objects/siteCovs2018.RData')
names(siteCovariate)

# stations info
stations <- read.csv("Stations2018.csv")
hunting <- read.csv("HuntingData2018.csv")
cameraRecords <- read.csv("RecordTable2018.csv")


############################################################################
########################### PICK A SPECIES #################################
####################### & FORMAT FOR UNMARKED ##############################
############################################################################

species <- peccary # pick which species to proceed with


# data as matrix
y <- as.matrix(species)
y <- (y[ order(as.numeric(row.names(y))), ]) #order matters


#################### IF YOU WANT TO POOL THE COLUMNS ###########################
# paste columns together, then if there's a 1 -> 1, if there's no 1 but a 0 -> 0, if only NAs -> NA
y <- as.matrix(species)
y <- (y[ order(as.numeric(row.names(y))), ]) #order matters

#### ***** CHANGE COMBINATION FACTOR HERE ***** ####
clumpEvery <- 10 
nClumpedColumns <- ncol(y)/clumpEvery
clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = 30) # bc there are 30 stations

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
############################################################################

# unmarked df for single species
ufo <- unmarkedFrameOccu(clumpedMatrix, # 'y' or 'combinedTime'
                        siteCovs = siteCovariate,
                        obsCovs = NULL)
peccaryUFO <- ufo
save(peccaryUFO, file = 'R Objects/PeccaryUFO.RData')

uncombinedUFO <- unmarkedFrameOccu(y, # 'y' or 'combinedTime'
                                   siteCovs = siteCovariate,
                                   obsCovs = NULL)

plot(ufo)
require(reshape2)
colnames(ufo@y) <- 1:ncol(ufo@y)
meltedComb <- melt(ufo@y)
meltedComb$value <- as.factor(meltedComb$value)
ggplot(meltedComb, aes(Var2, Var1, fill = value)) + 
  geom_tile(colour = "gray50") +
  scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
  scale_x_continuous(breaks = seq(0, ncol(ufo@y), by = 2)) +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  xlab("Time (1 unit = 8 days)") +
  ylab("Camera trap site") +
  theme_bw()

plot(uncombinedUFO)
# require(reshape2)
colnames(uncombinedUFO@y) <- 1:ncol(uncombinedUFO@y)
meltedUncomb <- melt(uncombinedUFO@y)
meltedUncomb$value <- as.factor(meltedUncomb$value)
ggplot(meltedUncomb, aes(Var2, Var1, fill = value)) + 
  geom_tile(colour = "gray50") +
  scale_fill_manual(values=c("gray75", "red3"), na.value="white", name = "") +
  scale_x_continuous(breaks = seq(0, ncol(uncombinedUFO@y), by = 4)) +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  xlab("Time (1 unit = 2 days)") +
  ylab("Camera trap site") +
  theme_bw()

############################################################################
############################## DETECTION ###################################
############################################################################

# vector of variables to consider for detection
match_variables <- c("Effort", "Habitat", "OnTrail")

# every possible combination of variables
combos <- sapply( seq(3), function(i) {
  as.list(as.data.frame(combn( x = match_variables, m = i)))
})
combos <- unlist(combos, recursive=FALSE)

# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))

detectionFormulas <- forms


############################################################################
############################## OCCUPANCY ###################################
############################################################################

# vector of variables 
match_variables <- c("Community", "River", "Habitat", 
                     "Hunting", "Trail.Distance", "HuntingIntensity")

# every possible combination of variables
combos <- sapply( seq(length(match_variables)), function(i) {
  as.list(as.data.frame(combn( x = match_variables, m = i)))
})
combos <- unlist(combos, recursive=FALSE)

# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))


occupancyFormulas <- forms


############################################################################
############################### MODELS #####################################
############################################################################

#### add in the null models
tempDF <- data.frame(detection = c("~1", 
                                   detectionFormulas,
                                   rep("~1", times = length(occupancyFormulas))),
                     occupancy = c("~1", 
                                   rep("~1", times = length(detectionFormulas)),
                                   occupancyFormulas))

##################### UNMARKED ########################
# every combination of occupancy and detection formulas
repeated <- data.frame(detection = rep(detectionFormulas, times = length(occupancyFormulas)),
                       occupancy = rep(occupancyFormulas, each = length(detectionFormulas)))
repeated <- rbind(repeated, tempDF)

allunmarkedFormulas <- paste(repeated$detection, repeated$occupancy, sep = " ")  

# run occupancy unmarked model for all those combinations
mods=list()
for(i in 1:length(allunmarkedFormulas)) {
  if (i%%20==0) print(paste("Model", i, "of", length(allunmarkedFormulas))) 
  test <- occu(formula(allunmarkedFormulas[[i]]), ufo)
  
  mods[[i]] <- occu(formula(allunmarkedFormulas[[i]]), ufo, 
                    control = 10000, # to get models to converge (https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html#common-errors)
                    starts = c(rep(0, length(test@opt$par)))) # setting all starting values to 0
}
# name the models
names(mods) <- 1:length(mods)
names(mods)[which(allunmarkedFormulas == "~1 ~1")] <- "Null"
mods[[which(allunmarkedFormulas == "~1 ~1")]]

peccaryModels <- mods
save(peccaryModels, file = 'R Objects/PeccaryAllModels.RData')

# remove all models with missing SE/z/p-value or that didn't converge
noMissingMods <- list()
for(i in 1:length(mods)){
  
  modSum <- summary(mods[[i]])
  
  if(anyNA(modSum$state$SE)==FALSE & mods[[i]]@opt$convergence != 1){
    noMissingMods[[i]] <- mods[[i]]
  } else {
    next
  }
}
noMissingMods[sapply(noMissingMods, is.null)] <- NULL

# test <- data.frame(Y = rep(NA, length(noMissingMods)))
# for(i in 1:length(noMissingMods)){
#   test[i,] <- ifelse(noMissingMods[[i]]@opt$convergence == 1, 1, 0)
# }
# sum(test)

# remove all models with ridiculous standard errors
noWackyMods <- list()
for(i in 1:length(noMissingMods)){
  
  modSum <- summary(noMissingMods[[i]])
  
  if((abs(modSum$state$SE[1]) < modSum$state$Estimate[1]) == TRUE){
    noWackyMods[[i]] <- noMissingMods[[i]]
  } else{
    next
  }
}

noWackyMods[sapply(noWackyMods, is.null)] <- NULL

peccaryNoWackyMods <- noWackyMods
save(peccaryNoWackyMods, file = 'R Objects/PeccaryNoWackyMods.RData')

# quick model selection
peccFitList <- fitList(noWackyMods)
peccaryModelSelection <- modSel(peccFitList)
peccaryModelSelection <- as(peccaryModelSelection, "data.frame")
save(peccaryModelSelection, file = 'R Objects/peccaryModelSelection.RData')


############################################################################
########################### MODEL SELECTION ################################
############################################################################

# Make a data frame to show the model names and their AICs
df <- data.frame(ModelName = NA,
                 AIC = NA)

for(i in 1:length(noWackyMods)){
  df[i,1]<- as.character(c(noWackyMods[[i]]@formula))
  df[i,2]<- noWackyMods[[i]]@AIC
}

unmarkedModels <- df[order(df$AIC),] # order by AIC

# calculate difference in AIC from #1 model
unmarkedModels$diffFromBest <- NA
for (i in 1:nrow(unmarkedModels)) {
  unmarkedModels$diffFromBest[i] <- unmarkedModels$AIC[i] - unmarkedModels$AIC[1]
}
#head(unmarkedModels, 10)

# take the best models
(topModels <- subset(unmarkedModels, diffFromBest <= 2)) # best occupancy models within 2 AIC of the lowest

# this is the best model
occu(formula(unmarkedModels[1,1]), ufo) 

# take these best models and put them into a list
bestMods <- list()
for(i in 1:nrow(topModels)) {
  test <- occu(formula(topModels[[i]]), ufo)
  
  bestMods[[i]] <- occu(formula(topModels[[i]]), ufo, 
                    control = 10000, # to get models to converge (https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html#common-errors)
                    starts = c(rep(0, length(test@opt$par)))) # setting all starting values to 0
  #bestMods[[i]] <- occu(formula(topModels[i,]$ModelName), ufo)
  print(bestMods[[i]])
} 

bestModsList <- bestMods
bestMods <- fitList(bestModsList)

allModelsOrdered <- list()
for(i in 1:nrow(unmarkedModels)) {
  test <- occu(formula(unmarkedModels[[i]]), ufo)
  
  allModelsOrdered[[i]] <- occu(formula(unmarkedModels[[i]]), ufo, 
                    control = 10000, # to get models to converge (https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html#common-errors)
                    starts = c(rep(0, length(test@opt$par)))) # setting all starting values to 0
  #allModelsOrdered[[i]] <- occu(formula(unmarkedModels[i,]$ModelName), ufo)
  #print(allModelsOrdered[[i]])
} 

# plot standard errors across the ranked models
orderedModelsEstimates <- data.frame(Estimate = rep(NA, times = length(allModelsOrdered)),
                                     SE = NA)
for (i in 1:length(allModelsOrdered)) {
  modSum <- summary(allModelsOrdered[[i]])
  orderedModelsEstimates$SE[i] <- modSum$state$SE[1]
  orderedModelsEstimates$Estimate[i] <- modSum$state$Estimate[1]
}

barplot(orderedModelsEstimates$SE/orderedModelsEstimates$Estimate, 
        xlab = "Models ordered by AIC (lowest -> highest)",
        ylab = "Standard error/estimate of model", 
        main = "Peccary")


############################################################################
##################### PLOTTING MODEL PREDICTIONS ###########################
############################################################################


############### COMMUNITY ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
df <- data.frame(Community = seq(min(siteCovariate$Community), 
                                 max(siteCovariate$Community), 
                                 length.out = 100),
                 Effort = mean(siteCovariate$Effort),
                 River = mean(siteCovariate$River),
                 Habitat = as.factor("Upland"), # most common habitat
                 Hunting = as.factor(1),
                 Trail.Distance = mean(siteCovariate$Trail.Distance),
                 OnTrail = as.factor(0), # more common to be off trail
                 HuntingIntensity = mean(siteCovariate$HuntingIntensity)
) 

# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)

# Put prediction, confidence interval, and covariate values together in a data frame
predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                  lower = unmarkedPred$lower,
                                  upper = unmarkedPred$upper,
                                  df)


# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = Community, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to community (scaled)", y = "Occupancy probability") +
  # ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 2.05)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
library(rphylopic)
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryComm <- predictionPlot + 
  add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))



############### RIVER ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
df <- data.frame(Community = mean(siteCovariate$Community),
                 Effort = mean(siteCovariate$Effort),
                 River = seq(min(siteCovariate$River), 
                             max(siteCovariate$River), 
                             length.out = 100),
                 Habitat = as.factor("Upland"), # most common habitat
                 Hunting = as.factor(1),
                 Trail.Distance = mean(siteCovariate$Trail.Distance),
                 OnTrail = as.factor(0) , # more common to be off trail
                 HuntingIntensity = mean(siteCovariate$HuntingIntensity)
) 

# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)

# Put prediction, confidence interval, and covariate values together in a data frame
predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                  lower = unmarkedPred$lower,
                                  upper = unmarkedPred$upper,
                                  df)

# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = River, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'blue2', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to river (scaled)", y = "Occupancy probability") +
  # ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 2.05)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
library(rphylopic)
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryRiver <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))



############### TRAIL DISTANCE ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
df <- data.frame(Community = mean(siteCovariate$Community),
                 Effort = mean(siteCovariate$Effort),
                 River = mean(siteCovariate$River),
                 Habitat = as.factor("Upland"), # most common habitat
                 Hunting = as.factor(1),
                 Trail.Distance = seq(min(siteCovariate$Trail.Distance), 
                                      max(siteCovariate$Trail.Distance), 
                                      length.out = 100),
                 OnTrail = as.factor(0) , # more common to be off trail
                 HuntingIntensity = mean(siteCovariate$HuntingIntensity)
) 

# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)

# Put prediction, confidence interval, and covariate values together in a data frame
predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                  lower = unmarkedPred$lower,
                                  upper = unmarkedPred$upper,
                                  df)

# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = Trail.Distance, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'green4', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to a trail (scaled)", y = "Occupancy probability") +
  #ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 2.05)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
library(rphylopic)
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryTrail <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))


############### HUNTING INTENSITY ##################
# First, set-up a new data frame to predict along a sequence of the covariate.
# Predicting requires all covariates, so hold the other covariates constant at their mean value
df <- data.frame(Community = mean(siteCovariate$Community),
                 Effort = mean(siteCovariate$Effort),
                 River = mean(siteCovariate$River),
                 Habitat = as.factor("Upland"), # most common habitat
                 Hunting = as.factor(1),
                 Trail.Distance = mean(siteCovariate$Trail.Distance),
                 OnTrail = as.factor(0), # more common to be off trail
                 HuntingIntensity = seq(min(siteCovariate$HuntingIntensity), 
                                        max(siteCovariate$HuntingIntensity), 
                                        length.out = 100)
) 

# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)

# Put prediction, confidence interval, and covariate values together in a data frame
predictionDataFrame <- data.frame(Predicted = unmarkedPred$Predicted,
                                  lower = unmarkedPred$lower,
                                  upper = unmarkedPred$upper,
                                  df)


# Plot the relationship
predictionPlot <- ggplot(predictionDataFrame, aes(x = HuntingIntensity, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'orange', alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Hunting intensity (scaled)", y = "Occupancy probability") +
  # ggtitle(label = "Collared peccary") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1), xlim = c(-1.85, 3)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, face = 'italic'))

# add an animal
library(rphylopic)
# Pecari tajacu
# Cervus elaphus
# Cuniculus paca

# Get a single image uuid for a species
uuid <- get_uuid(name = "Pecari tajacu", n = 1)
# Get the image for that uuid
img <- get_phylopic(uuid = uuid)

(peccaryHunting <- predictionPlot + 
    add_phylopic(img, alpha = 1, x = 1.3, y = 0.1, ysize = 0.5))










# figure panel
require(ggpubr)
arranged <- ggarrange(peccaryComm, peccaryRiver, peccaryTrail, peccaryHunting,
          #labels = c("A", "B", "C"),
          ncol = 4, nrow = 1)
annotate_figure(arranged, top = text_grob("Collared peccary", face = "bold", size = 16))


























###### average occupancy and detection across sites

df <- data.frame(Community = siteCovariate$Community,
                 Effort = siteCovariate$Effort,
                 River = siteCovariate$River,
                 Habitat = siteCovariate$Habitat,
                 Hunting = siteCovariate$Hunting,
                 Trail.Distance = siteCovariate$Trail.Distance,
                 OnTrail = siteCovariate$OnTrail) 

# Model-averaged prediction of occupancy and confidence interval
unmarkedPred <- predict(bestMods, type = 'state', new = df, appendData = TRUE)
mean(unmarkedPred$Predicted); mean(unmarkedPred$SE)
range(unmarkedPred$Predicted)

unmarkedPred <- predict(bestMods, type = 'det', new = df, appendData = TRUE)
mean(unmarkedPred$Predicted); mean(unmarkedPred$SE)
range(unmarkedPred$Predicted)
