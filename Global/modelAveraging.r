ManualModelAverage <- function(modelList, speciesName = "IDK") {
    # check if single model
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

        if (!singleModel && sumWeights > 0) {
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

    # create output dataframes
    # occupancy dataframe
    occupancyDF <- data.frame(
        Species = speciesName,
        PercentNaturalScaled = NA,
        PercentNaturalScaledSE = NA,
        PercentNaturalScaledUpper = NA,
        PercentNaturalScaledLower = NA,
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

    # detection dataframe
    detectionDF <- data.frame(
        Species = speciesName,
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

    # check if Community is a covariate
    hasCommunityInOcc <- any(grepl("^Community", rownames(occEstimates)))



    # fill occupancy dataframe with values
    for (param in rownames(occEstimates)) {
        # handle intercept as Zabalo if Community is a covariate
        if (param == "(Intercept)") {
            if (hasCommunityInOcc) {
                occupancyDF$CommunityZabalo <- occEstimates[param, "Estimate"]
                occupancyDF$CommunityZabaloSE <- occEstimates[param, "SE"]
                occupancyDF$CommunityZabaloLower <- occEstimates[param, "Lower"]
                occupancyDF$CommunityZabaloUpper <- occEstimates[param, "Upper"]
            }
            next
        }

        # match parameter names to dataframe columns
        if (param == "PercentNaturalScaled" && !is.na(occEstimates[param, "Estimate"])) {
            occupancyDF$PercentNaturalScaled <- occEstimates[param, "Estimate"]
            occupancyDF$PercentNaturalScaledSE <- occEstimates[param, "SE"]
            occupancyDF$PercentNaturalScaledLower <- occEstimates[param, "Lower"]
            occupancyDF$PercentNaturalScaledUpper <- occEstimates[param, "Upper"]
        } else if (param == "RainfallScaled" && !is.na(occEstimates[param, "Estimate"])) {
            occupancyDF$RainfallScaled <- occEstimates[param, "Estimate"]
            occupancyDF$RainfallScaledSE <- occEstimates[param, "SE"]
            occupancyDF$RainfallScaledLower <- occEstimates[param, "Lower"]
            occupancyDF$RainfallScaledUpper <- occEstimates[param, "Upper"]
        } else if (param == "DistToWater" && !is.na(occEstimates[param, "Estimate"])) {
            occupancyDF$DistToWater <- occEstimates[param, "Estimate"]
            occupancyDF$DistToWaterSE <- occEstimates[param, "SE"]
            occupancyDF$DistToWaterLower <- occEstimates[param, "Lower"]
            occupancyDF$DistToWaterUpper <- occEstimates[param, "Upper"]
        } else if (param == "TemperatureScaled" && !is.na(occEstimates[param, "Estimate"])) {
            occupancyDF$TemperatureScaled <- occEstimates[param, "Estimate"]
            occupancyDF$TemperatureScaledSE <- occEstimates[param, "SE"]
            occupancyDF$TemperatureScaledLower <- occEstimates[param, "Lower"]
            occupancyDF$TemperatureScaledUpper <- occEstimates[param, "Upper"]
        } else if (param == "DistToComm" && !is.na(occEstimates[param, "Estimate"])) {
            occupancyDF$DistToComm <- occEstimates[param, "Estimate"]
            occupancyDF$DistToCommSE <- occEstimates[param, "SE"]
            occupancyDF$DistToCommLower <- occEstimates[param, "Lower"]
            occupancyDF$DistToCommUpper <- occEstimates[param, "Upper"]
        } else if (grepl("^Community", param)) {
            # remove "Community" prefix, handle special cases like "San Pablo"
            commName <- gsub("^Community", "", param)

            # map the parameter name to the correct column name
            if (commName == "Remolino") {
                occupancyDF$CommunityRemolino <- occEstimates[param, "Estimate"]
                occupancyDF$CommunityRemolinoSE <- occEstimates[param, "SE"]
                occupancyDF$CommunityRemolinoLower <- occEstimates[param, "Lower"]
                occupancyDF$CommunityRemolinoUpper <- occEstimates[param, "Upper"]
            } else if (commName == "Sinangoe") {
                occupancyDF$CommunitySinangoe <- occEstimates[param, "Estimate"]
                occupancyDF$CommunitySinangoeSE <- occEstimates[param, "SE"]
                occupancyDF$CommunitySinangoeLower <- occEstimates[param, "Lower"]
                occupancyDF$CommunitySinangoeUpper <- occEstimates[param, "Upper"]
            } else if (commName == "San Pablo" || commName == "SanPablo") {
                occupancyDF$CommunitySanPablo <- occEstimates[param, "Estimate"]
                occupancyDF$CommunitySanPabloSE <- occEstimates[param, "SE"]
                occupancyDF$CommunitySanPabloLower <- occEstimates[param, "Lower"]
                occupancyDF$CommunitySanPabloUpper <- occEstimates[param, "Upper"]
            } else if (commName == "Siona") {
                occupancyDF$CommunitySiona <- occEstimates[param, "Estimate"]
                occupancyDF$CommunitySionaSE <- occEstimates[param, "SE"]
                occupancyDF$CommunitySionaLower <- occEstimates[param, "Lower"]
                occupancyDF$CommunitySionaUpper <- occEstimates[param, "Upper"]
            }
        }
    }

    # check if Community is a covariate
    hasCommunityInDet <- any(grepl("^Community", rownames(detEstimates)))

    # fill detection dataframe
    for (param in rownames(detEstimates)) {
        # handle intercept as Zabalo if Community is a covariate
        if (param == "(Intercept)") {
            if (hasCommunityInDet) {
                detectionDF$CommunityZabalo <- detEstimates[param, "Estimate"]
                detectionDF$CommunityZabaloSE <- detEstimates[param, "SE"]
                detectionDF$CommunityZabaloLower <- detEstimates[param, "Lower"]
                detectionDF$CommunityZabaloUpper <- detEstimates[param, "Upper"]
            }
            next
        }

        # match parameter names to dataframe columns
        if (param == "DaysEffortScaled") {
            detectionDF$DaysEffortScaled <- detEstimates[param, "Estimate"]
            detectionDF$DaysEffortScaledSE <- detEstimates[param, "SE"]
            detectionDF$DaysEffortScaledLower <- detEstimates[param, "Lower"]
            detectionDF$DaysEffortScaledUpper <- detEstimates[param, "Upper"]
        } else if (grepl("^Community", param)) {
            # remove "Community", handle special cases like "San Pablo"
            commName <- gsub("^Community", "", param)

            # map the parameter name to the correct column name
            if (commName == "Remolino") {
                detectionDF$CommunityRemolino <- detEstimates[param, "Estimate"]
                detectionDF$CommunityRemolinoSE <- detEstimates[param, "SE"]
                detectionDF$CommunityRemolinoLower <- detEstimates[param, "Lower"]
                detectionDF$CommunityRemolinoUpper <- detEstimates[param, "Upper"]
            } else if (commName == "Sinangoe") {
                detectionDF$CommunitySinangoe <- detEstimates[param, "Estimate"]
                detectionDF$CommunitySinangoeSE <- detEstimates[param, "SE"]
                detectionDF$CommunitySinangoeLower <- detEstimates[param, "Lower"]
                detectionDF$CommunitySinangoeUpper <- detEstimates[param, "Upper"]
            } else if (commName == "San Pablo" || commName == "SanPablo") {
                detectionDF$CommunitySanPablo <- detEstimates[param, "Estimate"]
                detectionDF$CommunitySanPabloSE <- detEstimates[param, "SE"]
                detectionDF$CommunitySanPabloLower <- detEstimates[param, "Lower"]
                detectionDF$CommunitySanPabloUpper <- detEstimates[param, "Upper"]
            } else if (commName == "Siona") {
                detectionDF$CommunitySiona <- detEstimates[param, "Estimate"]
                detectionDF$CommunitySionaSE <- detEstimates[param, "SE"]
                detectionDF$CommunitySionaLower <- detEstimates[param, "Lower"]
                detectionDF$CommunitySionaUpper <- detEstimates[param, "Upper"]
            }
        }
    }

    # return list with both dataframes and extra info
    result <- list(
        occupancy = occupancyDF,
        detection = detectionDF,
        # aicWeights = if(!singleModel) aicWeights else NULL,
        # deltaAIC = if(!singleModel) deltaAIC else NULL,
        unformattedOccupancy = occEstimates,
        unformattedDetection = detEstimates
    )

    return(result)
}
