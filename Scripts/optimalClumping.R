############################################################################
############################# LOAD DATA ####################################
############################################################################
# setwd("~/Documents/amazon/Data")
require(unmarked)

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

# function that takes the species occupancy data and tells you what the best clumping factor is
best_clumping_factor <- function(occupancyData){
  
  # make sure the first row and column are not row or column names
  if(sum(occupancyData[,1], na.rm = TRUE) > ncol(occupancyData) | # if the first row is column names
     sum(occupancyData[1,], na.rm = TRUE) > nrow(occupancyData)){ # if the first column in row names
    paste("Make sure your first row and first column are not row names or column names!")
  } else {
  
  y <- occupancyData
  allClumpedMatrices <- list()
  ncol <- ncol(occupancyData)
  
  for (i in 1:(ncol)) { # for every possible clumping factor i, make a clumped matrix
    clumpingFactor <- i
    nClumpedColumns <- ncol/clumpingFactor
    clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = 30) # bc there are 30 stations
    
    clumpStart <- seq(1, ncol, by = clumpingFactor) # the first column in the clump
    clumpEnd <- seq(clumpingFactor, ncol, by = clumpingFactor) # the last column in the clump
    
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
    
    allClumpedMatrices[[i]] <- clumpedMatrix
    
  } # all possible clumped matrices are made and stored in a list
  
  ### run model for each matrix
  allClumpedModels <- list()
  for (m in 1:length(allClumpedMatrices)) { 
    
    umf <- unmarkedFrameOccu(y = allClumpedMatrices[[m]])
    model <- occu(formula = ~1 ~1, data=umf)
    allClumpedModels[[m]] <- model # store all models

  }
  
  ### get SE for each possible model
  allSE <- data.frame(clumpingFactor = 1:ncol,
                      modelSE = NA)
  for (n in 1:length(allClumpedModels)) {
    allSE$modelSE[n] <- summary(allClumpedModels[[n]])$state$SE[1] # occupancy intercept SE
  }
  
  print(allSE)
  paste("The best clumping factor for this species is", allSE$clumpingFactor[which.min(allSE$modelSE)], ":)")
  
  } # end of if/else
}


############################################################################
########################### PICK A SPECIES #################################
####################### & FORMAT FOR UNMARKED ##############################
############################################################################

species <- paca # pick which species to proceed with
y <- as.matrix(species)
y <- (y[ order(as.numeric(row.names(y))), ]) #order matters

best_clumping_factor(y)
# peccary = 10
# brocket = 18
# paca = 1

