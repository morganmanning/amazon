ManualModelAverage <- function(modelList, speciesName = "IDK") {
    # check if only one model
    if (length(modelList) == 1) {
        singleModel <- TRUE
    } else {
        singleModel <- FALSE
    }
    nModels <- length(modelList)
    # extract AIC values
    aicValues <- numeric(nModels)
    for (i in 1:nModels) {
        aicValues[i] <- modelList[[i]]@AIC
    }
    # calculate AIC weights
    deltaAIC <- aicValues - min(aicValues)
    aicWeights <- exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))
    # get all unique parameter names from all models
    allOccupancyParams <- c()
    allDetectionParams <- c()
    for (i in 1:nModels) {
        occParams <- names(modelList[[i]]@estimates@estimates$state@estimates)
        detParams <- names(modelList[[i]]@estimates@estimates$det@estimates)
        allOccupancyParams <- unique(c(allOccupancyParams, occParams))
        allDetectionParams <- unique(c(allDetectionParams, detParams))
    }
    # model-averaged estimates
    occEstimates <- matrix(NA, nrow = length(allOccupancyParams), ncol = 4)
    rownames(occEstimates) <- allOccupancyParams
    colnames(occEstimates) <- c("Estimate", "SE", "Lower", "Upper")
    detEstimates <- matrix(NA, nrow = length(allDetectionParams), ncol = 4)
    rownames(detEstimates) <- allDetectionParams
    colnames(detEstimates) <- c("Estimate", "SE", "Lower", "Upper")
    # calculate z-value for confidence intervals
    zVal <- qnorm(1 - (1 - 0.95) / 2)
    # process occupancy parameters
    for (param in allOccupancyParams) {
        weightedEst <- 0
        weightedVar <- 0
        sumWeights <- 0
        for (i in 1:nModels) {
            # check if this parameter exists in this model
            paramNames <- names(modelList[[i]]@estimates@estimates$state@estimates)
            if (param %in% paramNames) {
                idx <- which(paramNames == param)
                estimate <- modelList[[i]]@estimates@estimates$state@estimates[idx]
                se <- sqrt(modelList[[i]]@estimates@estimates$state@covMat[idx, idx])
                if (singleModel) {
                    # for single model, just use the estimates
                    occEstimates[param, "Estimate"] <- estimate
                    occEstimates[param, "SE"] <- se
                } else {
                    # model averaging
                    weightedEst <- weightedEst + aicWeights[i] * estimate
                    sumWeights <- sumWeights + aicWeights[i]
                }
            }
        }
        if (!singleModel && sumWeights > 0) { # if not single model and parameter exists in at least one model
            # model averaging
            occEstimates[param, "Estimate"] <- weightedEst
            # calculate SE
            for (i in 1:nModels) {
                paramNames <- names(modelList[[i]]@estimates@estimates$state@estimates)
                if (param %in% paramNames) {
                    idx <- which(paramNames == param)
                    estimate <- modelList[[i]]@estimates@estimates$state@estimates[idx]
                    se <- sqrt(modelList[[i]]@estimates@estimates$state@covMat[idx, idx])
                    weightedVar <- weightedVar + aicWeights[i] * (se^2 + (estimate - weightedEst)^2)
                }
            }
            occEstimates[param, "SE"] <- sqrt(weightedVar)
        }
        # calculate confidence intervals
        if (!is.na(occEstimates[param, "Estimate"])) {
            occEstimates[param, "Lower"] <- occEstimates[param, "Estimate"] - zVal * occEstimates[param, "SE"]
            occEstimates[param, "Upper"] <- occEstimates[param, "Estimate"] + zVal * occEstimates[param, "SE"]
        }
    }
    # process detection parameters
    for (param in allDetectionParams) {
        weightedEst <- 0
        weightedVar <- 0
        sumWeights <- 0
        for (i in 1:nModels) {
            # check if this parameter exists in this model
            paramNames <- names(modelList[[i]]@estimates@estimates$det@estimates)
            if (param %in% paramNames) {
                idx <- which(paramNames == param)
                estimate <- modelList[[i]]@estimates@estimates$det@estimates[idx]
                se <- sqrt(modelList[[i]]@estimates@estimates$det@covMat[idx, idx])
                if (singleModel) {
                    # for single model, just use the estimates
                    detEstimates[param, "Estimate"] <- estimate
                    detEstimates[param, "SE"] <- se
                } else {
                    # model averaging
                    weightedEst <- weightedEst + aicWeights[i] * estimate
                    sumWeights <- sumWeights + aicWeights[i]
                }
            }
        }
        if (!singleModel && sumWeights > 0) {
            # model averaging
            detEstimates[param, "Estimate"] <- weightedEst
            # calculate SE
            for (i in 1:nModels) {
                paramNames <- names(modelList[[i]]@estimates@estimates$det@estimates)
                if (param %in% paramNames) {
                    idx <- which(paramNames == param)
                    estimate <- modelList[[i]]@estimates@estimates$det@estimates[idx]
                    se <- sqrt(modelList[[i]]@estimates@estimates$det@covMat[idx, idx])
                    weightedVar <- weightedVar + aicWeights[i] * (se^2 + (estimate - weightedEst)^2)
                }
            }
            detEstimates[param, "SE"] <- sqrt(weightedVar)
        }
        # calculate confidence intervals
        if (!is.na(detEstimates[param, "Estimate"])) {
            detEstimates[param, "Lower"] <- detEstimates[param, "Estimate"] - zVal * detEstimates[param, "SE"]
            detEstimates[param, "Upper"] <- detEstimates[param, "Estimate"] + zVal * detEstimates[param, "SE"]
        }
    }
        
    # check if Community is a covariate
    hasCommunityInOcc <- any(grepl("^Community", rownames(occEstimates)))
    hasCommunityInDet <- any(grepl("^Community", rownames(detEstimates)))
    
    # create occupancy dataframe dynamically
    occColNames <- c("Species")
    occColValues <- list(Species = speciesName)
    
    for (param in rownames(occEstimates)) {
        # handle intercept as Zabalo if Community is a covariate
        if (param == "(Intercept)") {
            if (hasCommunityInOcc) {
                paramName <- "CommunityZabalo"
            } else {
                next  # skip intercept if no Community covariate
            }
        } else {
            paramName <- param
        }
        
        # add columns for this parameter
        occColNames <- c(occColNames, paramName, 
                        paste0(paramName, "SE"),
                        paste0(paramName, "Upper"),
                        paste0(paramName, "Lower"))
        
        occColValues[[paramName]] <- occEstimates[param, "Estimate"]
        occColValues[[paste0(paramName, "SE")]] <- occEstimates[param, "SE"]
        occColValues[[paste0(paramName, "Upper")]] <- occEstimates[param, "Upper"]
        occColValues[[paste0(paramName, "Lower")]] <- occEstimates[param, "Lower"]
    }
    
    occupancyDF <- as.data.frame(occColValues, stringsAsFactors = FALSE)
    
    # create detection dataframe dynamically
    detColNames <- c("Species")
    detColValues <- list(Species = speciesName)
    
    for (param in rownames(detEstimates)) {
        # handle intercept as Zabalo if Community is a covariate
        if (param == "(Intercept)") {
            if (hasCommunityInDet) {
                paramName <- "CommunityZabalo"
            } else {
                next  # skip intercept if no Community covariate
            }
        } else {
            paramName <- param
        }
        
        # add columns for this parameter
        detColNames <- c(detColNames, paramName, 
                        paste0(paramName, "SE"),
                        paste0(paramName, "Upper"),
                        paste0(paramName, "Lower"))
        
        detColValues[[paramName]] <- detEstimates[param, "Estimate"]
        detColValues[[paste0(paramName, "SE")]] <- detEstimates[param, "SE"]
        detColValues[[paste0(paramName, "Upper")]] <- detEstimates[param, "Upper"]
        detColValues[[paste0(paramName, "Lower")]] <- detEstimates[param, "Lower"]
    }
    
    detectionDF <- as.data.frame(detColValues, stringsAsFactors = FALSE)
    
    # return list with both dataframes and extra info
    result <- list(
        occupancy = occupancyDF,
        detection = detectionDF,
        unformattedOccupancy = occEstimates,
        unformattedDetection = detEstimates
    )
    return(result)
}