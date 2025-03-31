############################################################################
############################# LOAD DATA ####################################
############################################################################
# setwd("~/Documents/amazon/Data")
require(unmarked)

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
    clumpedMatrix <- matrix(0, ncol = nClumpedColumns, nrow = nrow(occupancyData))
    
    clumpStart <- seq(1, ncol, by = clumpingFactor) # the first column in the clump
    clumpEnd <- seq(clumpingFactor, ncol, by = clumpingFactor) # the last column in the clump
    
    ### make the clumped matrix
    for (k in 1:nrow(clumpedMatrix)){
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
  for (m in 1:(length(allClumpedMatrices))) { 
    
    umf <- unmarkedFrameOccu(y = allClumpedMatrices[[m]])
    model <- occu(formula = ~1 ~1, data = umf)
    allClumpedModels[[m]] <- model # store all models

  }
  
  ### get SE for each possible model
  allSE <- data.frame(clumpingFactor = 1:ncol,
                      modelSE = NA)
  for (n in 1:length(allClumpedModels)) {
    # if there is an error when trying to extract the SE, then set it to NA
    # tryCatch is a function that allows you to catch errors in R
    # it will run the code in the first argument, and if there is an error, it will run the code in the second argument
    # in this case, if there is an error, it will set the SE to NA
    result <- tryCatch(
        {
            summary(allClumpedModels[[n]])$state$SE[1] # occupancy intercept SE
        },
      error = function(e) e
    )
    if (inherits(result, "error")) {
        allSE$modelSE[n] <- NA
    } else {
        allSE$modelSE[n] <- result
    }
    
    # allSE$modelSE[n] <- summary(allClumpedModels[[n]])$state$SE[1] # occupancy intercept SE
  }
  
  print(allSE)
  # paste("The best clumping factor for this species is", allSE$clumpingFactor[which.min(allSE$modelSE)], ":)")
  return(allSE$clumpingFactor[which.min(allSE$modelSE)])
  
  } # end of if/else
}

