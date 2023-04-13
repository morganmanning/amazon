###### PARAMETER ESTIMATION WORK SHOP ######

# setwd("~/Library/CloudStorage/Dropbox/UF/Spring 2023/WIS 6934 Parameter Estimation/Data")

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
siteCovariate <- read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$Community <- scale(siteCovariate$Community/1000)
siteCovariate$River <- scale(siteCovariate$River/1000)
siteCovariate$Effort <- scale(siteCovariate$Effort)
siteCovariate$OnTrail <- as.factor(ifelse(siteCovariate$Trail.Distance == 0, 1, 0))
siteCovariate$Trail.Distance <- scale(siteCovariate$Trail.Distance)
siteCovariate$Station <- NULL
siteCovariate$RR <- NULL
siteCovariate$CR <- NULL


# stations info
stations <- read.csv("Stations2018.csv")



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
clumpEvery <- 4

out <- sapply(seq(1, ncol(y), by = clumpEvery), function(i)  # every 4 time steps are clumped
  do.call(paste0, as.data.frame(y[, i:(i+(clumpEvery-1))]))) 


combinedTime <- matrix(nrow = nrow(out), ncol = ncol(out))
for(i in 1:ncol(out)){
  for(j in 1:nrow(out)){
    combinedTime[j,i] <- ifelse(grepl('1', out[j,i], fixed = TRUE)==TRUE, 1, 
                                ifelse(grepl('0', out[j,i], fixed = TRUE)==TRUE, 0, NA))
  }
}

#combinedTime <- combinedTime[-c(6,27),] # removing sites with lots of NAs
#siteCovariate <- siteCovariate[-c(6,27),]
############################################################################

# unmarked df for single species
ufo <- unmarkedFrameOccu(combinedTime, # 'y' or 'combinedTime'
                        siteCovs = siteCovariate,
                        obsCovs = NULL)

plot(ufo)


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
match_variables <- c("Community", "River", "Habitat", "Hunting", "Trail.Distance")

# every possible combination of variables
combos <- sapply( seq(5), function(i) {
  as.list(as.data.frame(combn( x = match_variables, m = i)))
})
combos <- unlist(combos, recursive=FALSE)

# all combinations of variables into formulas
forms <- sapply(combos, function(x) paste("~ ", paste(x, collapse="+"), sep = ""))


occupancyFormulas <- forms


############################################################################
############################### MODELS #####################################
############################################################################


##################### UNMARKED ########################
# every combination of occupancy and detection formulas
repeated <- data.frame(detection = rep(detectionFormulas, times = length(occupancyFormulas)),
                       occupancy = rep(occupancyFormulas, each = length(detectionFormulas)))
allunmarkedFormulas <- paste(repeated$detection, repeated$occupancy, sep = " ")  

# run occupancy unmarked model for all those combinations
mods=list()
for(i in 1:length(allunmarkedFormulas)) {
  mods[[i]] <- occu(formula(allunmarkedFormulas[[i]]), ufo)
}

# a lot of models are wonky
mods[[10]] # mods[[10]] for clumping times by 4; mods[[13]] for clumping times by 8
occu(~1 ~1, ufo) # but null model runs fine


# remove all models with missing SE/z/p-value
noMissingMods <- list()
for(i in 1:length(mods)){
  
  modSum <- summary(mods[[i]])
  
  if(anyNA(modSum$state$SE)==FALSE){
    noMissingMods[[i]] <- mods[[i]]
  } else {
    next
  }
}

noMissingMods[sapply(noMissingMods, is.null)] <- NULL


############################################################################
########################### MODEL SELECTION ################################
############################################################################

# Make a data frame to show the model names and their AICs
df <- data.frame(ModelName = NA,
                 AIC = NA)

for(i in 1:length(noMissingMods)){
  df[i,1]<- as.character(c(noMissingMods[[i]]@formula))
  df[i,2]<- noMissingMods[[i]]@AIC
}

unmarkedModels <- df[order(df$AIC),] # order by AIC

# calculate difference in AIC from #1 model
unmarkedModels$diffFromBest <- NA
for (i in 1:nrow(unmarkedModels)) {
  unmarkedModels$diffFromBest[i] <- unmarkedModels$AIC[i] - unmarkedModels$AIC[1]
}
head(unmarkedModels)

# take the best models
topModels <- subset(unmarkedModels, diffFromBest <= 2) # best occupancy models within 2 AIC of the lowest

# this is the best model
occu(formula(unmarkedModels[1,1]), ufo) # doesn't converge

# take these best models and put them into a list
bestMods <- list()
for(i in 1:nrow(topModels)) {
  bestMods[[i]] <- occu(formula(unmarkedModels[i,]$ModelName), ufo)
  print(bestMods[[i]])
} 
# 7 of the 21 best models didn't converge when time steps were reduced to 12
# 5 of the 6 best models didn't converge when time steps were reduced to 6

bestModsList <- bestMods
bestMods <- fitList(bestModsList)



############################################################################
##################### PLOTTING MODEL PREDICTIONS ###########################
############################################################################

require(AICcmodavg)

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
                 OnTrail = as.factor(0) # more common to be off trail
) 

# Model-averaged prediction of occupancy and confidence interval
occupancyModelAverage <- modavgPred(bestModsList, 
                                    parm.type = "psi", # psi = occupancy
                                    newdata = df)[c("mod.avg.pred",
                                                    "lower.CL",
                                                    "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
predictionDataFrame <- data.frame(Predicted = occupancyModelAverage$mod.avg.pred,
                                  lower = occupancyModelAverage$lower.CL,
                                  upper = occupancyModelAverage$upper.CL,
                                  df)


# Plot the relationship
(predictionPlot <- ggplot(predictionDataFrame, aes(x = Community, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(linewidth = 1) +
  labs(x = "Distance to community (km, scaled)", y = "Occupancy probability") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black")))




