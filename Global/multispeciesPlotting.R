##### PLOTTING MULTISPECIES MODELS #####
setwd("/Users/morganmanning/Documents/amazon/Global/Data")
setwd("~/Documents/amazon/Global/Data")

################################################################################
# ------------------------------ START UP -------------------------------------#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# load necessary packages
require(dplyr)
require(lubridate)
require(unmarked)
require(ggplot2)
require(rphylopic)
require(knitr)
require(kableExtra)
require(stringr)
require(gridExtra)
require(ggpubr)
require(reshape2)


# load the data
load("R Objects/multispeciesModels.RData")


################################################################################
###################### MARGINAL OCCUPANCY PREDICTIONS ##########################
################################################################################
# make a dataframe with predictions for all species
all_predictions <- data.frame()

for (i in 1:length(casualNames)) {
    # get predictions for each species
    preds <- predict(null_multispecies_model, type = "state", species = casualNames[i])
    all_predictions <- rbind(all_predictions, preds[1, ])
}
all_predictions$Species <- casualNames

# plot null model predictions for occupancy
plot(1:length(species), all_predictions$Predicted,
    # ylim = c(0.1, 0.4),
    xlim = c(0.5, 2.5),
    pch = 19, cex = 1.5, xaxt = "n",
    xlab = "", ylab = "Marginal occupancy and 95% CI",
    main = "Null model predictions for occupancy",
    ylim = c(0, 1)
)
axis(1, at = 1:length(species), labels = all_predictions$Species)

# CIs
top <- 0.1
for (i in 1:length(species)) {
    segments(i, all_predictions$lower[i], i, all_predictions$upper[i])
    segments(i - top, all_predictions$lower[i], i + top)
    segments(i - top, all_predictions$upper[i], i + top)
}

# ggplot version
all_predictions$Species <- factor(all_predictions$Species, levels = casualNames) # so plotting doesn't alphabetize species
ggplot(all_predictions, aes(x = Species, y = Predicted)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    ylim(0, 1) +
    labs(
        x = "Species",
        y = "Marginal occupancy and 95% CI"
    ) +
    theme_bw()
ggsave("../Figures/MultispeciesModeling/NullModelPredictions.png",
    width = 5, height = 5
)








################################################################################
#################### CONDITIONAL OCCUPANCY PREDICTIONS #########################
################################################################################

# make a dataframe with column of species 1 and column of species 2 with all combinations of species
# conditional_occupancy <- data.frame(
#     Species1 = combn(casualNames, 2)[1,],
#     Species2 = combn(casualNames, 2)[2, ],
#     Present1_Present2 = NA,
#     Present1_Absent2 = NA
# )

conditional_occupancy <- data.frame(
    Species1 = expand.grid(casualNames, casualNames)$Var2,
    Species2 = expand.grid(casualNames, casualNames)$Var1,
    Present1_Present2 = NA,
    PresentSE = NA,
    PresentLower = NA,
    PresentUpper = NA,
    Present1_Absent2 = NA,
    AbsentSE = NA,
    AbsentLower = NA,
    AbsentUpper = NA
)
conditional_occupancy <- conditional_occupancy[conditional_occupancy$Species1 != conditional_occupancy$Species2, ]

natural_conditional_occupancy <- data.frame(
    Species1 = expand.grid(casualNames, casualNames)$Var2,
    Species2 = expand.grid(casualNames, casualNames)$Var1,
    Present1_Present2 = NA,
    PresentSE = NA,
    PresentLower = NA,
    PresentUpper = NA,
    Present1_Absent2 = NA,
    AbsentSE = NA,
    AbsentLower = NA,
    AbsentUpper = NA
)
natural_conditional_occupancy <- natural_conditional_occupancy[natural_conditional_occupancy$Species1 != natural_conditional_occupancy$Species2, ]

best_conditional_occupancy <- data.frame(
    Species1 = expand.grid(casualNames, casualNames)$Var2,
    Species2 = expand.grid(casualNames, casualNames)$Var1,
    Present1_Present2 = NA,
    PresentSE = NA,
    PresentLower = NA,
    PresentUpper = NA,
    Present1_Absent2 = NA,
    AbsentSE = NA,
    AbsentLower = NA,
    AbsentUpper = NA
)
best_conditional_occupancy <- best_conditional_occupancy[best_conditional_occupancy$Species1 != best_conditional_occupancy$Species2, ]


for (i in 1:nrow(conditional_occupancy)) {
    # get predictions for each row
    Species1 <- conditional_occupancy$Species1[i]
    Species2 <- conditional_occupancy$Species2[i]

    # with both species present
    conditional_occupancy[i, c(3:6)] <- predict(null_multispecies_model,
        type = "state",
        species = Species1,
        cond = Species2
    )[1, ] # pull just the first row since rows are duplicates
    natural_conditional_occupancy[i, c(3:6)] <- colMeans(predict(natural_mod_penalty,
        type = "state",
        species = Species1,
        cond = Species2
    ))
    best_conditional_occupancy[i, c(3:6)] <- colMeans(predict(best_mod_penalty,
        type = "state",
        species = Species1,
        cond = Species2
    ))

    # with species 2 absent
    conditional_occupancy[i, c(7:10)] <- predict(null_multispecies_model,
        type = "state",
        species = Species1,
        cond = paste0("-", Species2)
    )[1, ]
    natural_conditional_occupancy[i, c(7:10)] <- colMeans(predict(natural_mod_penalty,
        type = "state",
        species = Species1,
        cond = paste0("-", Species2)
    ))
    best_conditional_occupancy[i, c(7:10)] <- colMeans(predict(best_mod_penalty,
        type = "state",
        species = Species1,
        cond = paste0("-", Species2)
    ))
}
### NULL MODEL
# make a ggplot panel figure for each species with the conditional occupancy with each other species
for (i in 1:length(casualNames)) {
    # get the species
    speciesInQuestion <- casualNames[i]

    # get the conditional occupancy for that species
    species_conditional_occupancy <- conditional_occupancy[conditional_occupancy$Species1 == speciesInQuestion, ]

    # plot each interaction
    species_plot_list <- list()
    for (j in 1:nrow(species_conditional_occupancy)) {
        plot_df <- data.frame(
            Status = c("Present", "Absent"),
            Predicted = c(species_conditional_occupancy$Present1_Present2[j], species_conditional_occupancy$Present1_Absent2[j]),
            SE = c(species_conditional_occupancy$PresentSE[j], species_conditional_occupancy$AbsentSE[j]),
            lower = c(species_conditional_occupancy$PresentLower[j], species_conditional_occupancy$AbsentLower[j]),
            upper = c(species_conditional_occupancy$PresentUpper[j], species_conditional_occupancy$AbsentUpper[j])
        )
        # plot plot_df
        species_plot_list[[j]] <- ggplot(plot_df, aes(x = Status, y = Predicted)) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
            ylim(0, 1) +
            labs(
                x = paste0(str_to_title(species_conditional_occupancy$Species2[j]), " status"),
                y = paste0(str_to_title(speciesInQuestion), " occupancy and 95% CI")
            ) +
            theme_bw()
    }

    # make the plot
    n <- length(species_plot_list)
    nCol <- floor(sqrt(n))
    p <- do.call("grid.arrange", c(species_plot_list, ncol = nCol))
    annotate_figure(p, top = text_grob(paste0(commonNames[i], " Interactions (Null model)"),
        face = "bold", size = 14
    ))
    ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[i], "_NullInteractions.png"),
        width = 10, height = 8
    )
}


#### NATURAL PENALIZED MODEL
# make a ggplot panel figure for each species with the conditional occupancy with each other species
for (i in 1:length(casualNames)) {
    # get the species
    speciesInQuestion <- casualNames[i]

    # get the conditional occupancy for that species
    species_conditional_occupancy <- natural_conditional_occupancy[natural_conditional_occupancy$Species1 == speciesInQuestion, ]

    # plot each interaction
    species_plot_list <- list()
    for (j in 1:nrow(species_conditional_occupancy)) {
        plot_df <- data.frame(
            Status = c("Present", "Absent"),
            Predicted = c(species_conditional_occupancy$Present1_Present2[j], species_conditional_occupancy$Present1_Absent2[j]),
            SE = c(species_conditional_occupancy$PresentSE[j], species_conditional_occupancy$AbsentSE[j]),
            lower = c(species_conditional_occupancy$PresentLower[j], species_conditional_occupancy$AbsentLower[j]),
            upper = c(species_conditional_occupancy$PresentUpper[j], species_conditional_occupancy$AbsentUpper[j])
        )
        # plot plot_df
        species_plot_list[[j]] <- ggplot(plot_df, aes(x = Status, y = Predicted)) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
            ylim(0, 1) +
            labs(
                x = paste0(str_to_title(species_conditional_occupancy$Species2[j]), " status"),
                y = paste0(str_to_title(speciesInQuestion), " occupancy and 95% CI")
            ) +
            theme_bw()
    }

    # make the plot
    n <- length(species_plot_list)
    nCol <- floor(sqrt(n))
    p <- do.call("grid.arrange", c(species_plot_list, ncol = nCol))
    annotate_figure(p, top = text_grob(paste0(commonNames[i], " Interactions (Penalized natural area model)"),
        face = "bold", size = 14
    ))
    ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[i], "_PenalizedNaturalInteractions.png"),
        width = 10, height = 8
    )
}


#### PENALIZED BEST MODEL
# make a ggplot panel figure for each species with the conditional occupancy with each other species
for (i in 1:length(casualNames)) {
    # get the species
    speciesInQuestion <- casualNames[i]

    # get the conditional occupancy for that species
    species_conditional_occupancy <- best_conditional_occupancy[best_conditional_occupancy$Species1 == speciesInQuestion, ]

    # plot each interaction
    species_plot_list <- list()
    for (j in 1:nrow(species_conditional_occupancy)) {
        plot_df <- data.frame(
            Status = c("Present", "Absent"),
            Predicted = c(species_conditional_occupancy$Present1_Present2[j], species_conditional_occupancy$Present1_Absent2[j]),
            SE = c(species_conditional_occupancy$PresentSE[j], species_conditional_occupancy$AbsentSE[j]),
            lower = c(species_conditional_occupancy$PresentLower[j], species_conditional_occupancy$AbsentLower[j]),
            upper = c(species_conditional_occupancy$PresentUpper[j], species_conditional_occupancy$AbsentUpper[j])
        )
        # plot plot_df
        species_plot_list[[j]] <- ggplot(plot_df, aes(x = Status, y = Predicted)) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
            ylim(0, 1) +
            labs(
                x = paste0(str_to_title(species_conditional_occupancy$Species2[j]), " status"),
                y = paste0(str_to_title(speciesInQuestion), " occupancy and 95% CI")
            ) +
            theme_bw()
    }

    # make the plot
    n <- length(species_plot_list)
    nCol <- floor(sqrt(n))
    p <- do.call("grid.arrange", c(species_plot_list, ncol = nCol))
    annotate_figure(p, top = text_grob(paste0(commonNames[i], " Interactions (Penalized best model)"),
        face = "bold", size = 14
    ))
    ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[i], "_PenalizedBestInteractions.png"),
        width = 10, height = 8
    )
}






################################################################################
########################### COVARIATE PREDICTIONS ##############################
################################################################################

covariates <- names(siteCovs(umf))
allCommunities <- unique(siteCovs(umf)$Community)
N <- 100
dfTemplate <- data.frame(
    X = rep(1, N),
    Station = rep("SGE1", N), # picked because it's kind of a middle community
    Community = rep("Sinangoe", N), # picked because it's kind of a middle community
    Rainfall = mean(siteCovs(umf)$Rainfall),
    RainfallScaled = mean(siteCovs(umf)$RainfallScaled),
    percentNatural = mean(siteCovs(umf)$percentNatural),
    DistToWater = mean(siteCovs(umf)$DistToWater),
    Temperature = mean(siteCovs(umf)$Temperature),
    TemperatureScaled = mean(siteCovs(umf)$TemperatureScaled),
    DistToComm = mean(siteCovs(umf)$DistToComm),
    NearestCommunity = rep("Sinangoe", N),
    Year = as.factor(rep(2022, N))
)

# extract the words from the best multispecies model state formulas
prediction_plots_list <- list() # plots per covariate included in the best models
covs <- unique(unlist(strsplit(gsub("1", "", gsub("[ ~]", "", best_mod_penalty@stateformulas)), "[ +~]")))
covs <- covs[covs != "Year" & covs != "Community"]
for (i in 1:length(covs)) {
    # make a prediction dataframe across the gradient of each covariate
    dfEdited <- dfTemplate
    covariateInQuestion <- covs[i]
    dfEdited[, covariateInQuestion] <- seq(min(siteCovs(umf)[, covariateInQuestion]),
        max(siteCovs(umf)[, covariateInQuestion]),
        length.out = N
    )
    plotting_df <- data.frame(
        Species = rep(NA, times = N),
        Covariate = NA,
        Predicted = NA,
        SE = NA,
        Lower = NA,
        Upper = NA
    )
    all_plotting_df <- list()
    for (j in 1:length(commonNames)) {
        # get the species
        speciesInQuestion <- casualNames[j]

        # get the predictions for the covariate
        preds <- predict(best_mod_penalty,
            type = "state",
            species = speciesInQuestion,
            newdata = dfEdited
        )
        plotting_df$Gradient <- dfEdited[, covariateInQuestion]
        plotting_df$Predicted <- preds$Predicted
        plotting_df$SE <- preds$SE
        plotting_df$Species <- speciesInQuestion
        plotting_df$Covariate <- covariateInQuestion
        plotting_df$Lower <- preds$lower
        plotting_df$Upper <- preds$upper
        all_plotting_df[[j]] <- plotting_df
    }
    big_plotting_df <- do.call("rbind", all_plotting_df)

    # plot predictions with ggplot with a line for each species
    big_plotting_df$Species <- factor(big_plotting_df$Species, levels = casualNames) # so plotting doesn't alphabetize species
    p <- ggplot(big_plotting_df, aes(x = Gradient, y = Predicted, color = Species)) +
        ylim(0, 1) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Species), alpha = 0.2, color = NA) +
        geom_line() +
        labs(
            x = str_to_title(covariateInQuestion),
            y = "Occupancy and 95% CI"
        ) +
        theme_bw()
    prediction_plots_list[[i]] <- p
    ggsave(
        file = paste0("../Figures/MultispeciesModeling/", covs[i], "_BestPredictions.png"),
        plot = p,
        width = 7, height = 7
    )
}
names(prediction_plots_list) <- covs


# extract the words from the best multispecies model state formulas
natural_prediction_plots_list <- list() # plots per covariate included in the best models
covs <- unique(unlist(strsplit(gsub("1", "", gsub("[ ~]", "", natural_mod_penalty@stateformulas)), "[ +~]")))
covs <- covs[covs != "Year" & covs != "Community"]
for (i in 1:length(covs)) {
    # make a prediction dataframe across the gradient of each covariate
    dfEdited <- dfTemplate
    covariateInQuestion <- covs[i]
    dfEdited[, covariateInQuestion] <- seq(min(siteCovs(umf)[, covariateInQuestion]),
        max(siteCovs(umf)[, covariateInQuestion]),
        length.out = N
    )
    plotting_df <- data.frame(
        Species = rep(NA, times = N),
        Covariate = NA,
        Predicted = NA,
        SE = NA,
        Lower = NA,
        Upper = NA
    )
    all_plotting_df <- list()
    for (j in 1:length(commonNames)) {
        # get the species
        speciesInQuestion <- casualNames[j]

        # get the predictions for the covariate
        preds <- predict(natural_mod_penalty,
            type = "state",
            species = speciesInQuestion,
            newdata = dfEdited
        )
        plotting_df$Gradient <- dfEdited[, covariateInQuestion]
        plotting_df$Predicted <- preds$Predicted
        plotting_df$SE <- preds$SE
        plotting_df$Species <- speciesInQuestion
        plotting_df$Covariate <- covariateInQuestion
        plotting_df$Lower <- preds$lower
        plotting_df$Upper <- preds$upper
        all_plotting_df[[j]] <- plotting_df
    }
    big_plotting_df <- do.call("rbind", all_plotting_df)

    # plot predictions with ggplot with a line for each species
    big_plotting_df$Species <- factor(big_plotting_df$Species, levels = casualNames) # so plotting doesn't alphabetize species
    p <- ggplot(big_plotting_df, aes(x = Gradient, y = Predicted, color = Species)) +
        ylim(0, 1) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Species), alpha = 0.2, color = NA) +
        geom_line() +
        labs(
            x = str_to_title(covariateInQuestion),
            y = "Occupancy and 95% CI"
        ) +
        theme_bw()
    natural_prediction_plots_list[[i]] <- p
    ggsave(
        file = paste0("../Figures/MultispeciesModeling/", covs[i], "_NaturalPredictions.png"),
        plot = p,
        width = 7, height = 7
    )
}
names(natural_prediction_plots_list) <- covs

# extract the words from the best multispecies model state formulas
community_prediction_plots_list <- list() # plots per covariate included in the best models
covs <- unique(unlist(strsplit(gsub("1", "", gsub("[ ~]", "", community_mod_penalty@stateformulas)), "[ +~]")))
for (i in 1:length(covs)) {
    # make a prediction dataframe across the gradient of each covariate
    dfEdited <- dfTemplate
    covariateInQuestion <- covs[i]
    orderedCommunities <- factor(unique(siteCovs(umf)[, covariateInQuestion]),
        levels = c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    )
    dfEdited[, covariateInQuestion] <- rep(orderedCommunities,
        each = N / length(unique(orderedCommunities))
    )

    plotting_df <- data.frame(
        Species = rep(NA, times = N),
        Covariate = NA,
        Predicted = NA,
        SE = NA,
        Lower = NA,
        Upper = NA
    )
    all_plotting_df <- list()
    for (j in 1:length(commonNames)) {
        # get the species
        speciesInQuestion <- casualNames[j]

        # get the predictions for the covariate
        preds <- predict(community_mod_penalty,
            type = "state",
            species = speciesInQuestion,
            newdata = dfEdited
        )
        plotting_df$Gradient <- dfEdited[, covariateInQuestion]
        plotting_df$Predicted <- preds$Predicted
        plotting_df$SE <- preds$SE
        plotting_df$Species <- speciesInQuestion
        plotting_df$Covariate <- covariateInQuestion
        plotting_df$Lower <- preds$lower
        plotting_df$Upper <- preds$upper
        all_plotting_df[[j]] <- plotting_df
    }
    big_plotting_df <- do.call("rbind", all_plotting_df)
    big_plotting_df <- big_plotting_df[!duplicated(big_plotting_df), ]
    big_plotting_df$Community <- gsub("Zabalo", "Zábalo", x = big_plotting_df$Gradient)
    big_plotting_df$Community <- factor(big_plotting_df$Community,
        levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    )


    # plot predictions with ggplot with a line for each species
    # make labels per covariate
    colors <- c(
        "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
        "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
    )

    # make axis titles per potential species
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

    # plot it
    dodge <- position_dodge(width = 0.3)
    big_plotting_df$Species <- factor(big_plotting_df$Species, levels = casualNames) # so plotting doesn't alphabetize species
    p <- ggplot(big_plotting_df, aes(
        x = Species,
        y = Predicted,
        color = Community
    )) +
        geom_point(aes(color = Community), position = dodge, size = 1.5) +
        geom_errorbar(
            aes(
                ymin = Predicted - SE,
                ymax = Predicted + SE,
                color = Community
            ),
            position = dodge, width = 0.15, linewidth = .5
        ) +
        # scale_color_manual(values = c("darkorange", "royalblue", "green3", "yellow3")) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        scale_x_discrete(labels = c(paca, agouti)) +
        labs(x = "Species", y = "Occupancy probability estimate (and SE)") +
        ylim(c(0, 1)) +
        theme_classic() +
        theme(
            text = element_text(family = "Times", colour = "black"),
            axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 45, vjust = 0.60),
            legend.title = element_blank(),
            legend.position = "top",
            axis.title.x = element_blank(),
            panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)
        ) +
        add_phylopic(pacaPic, alpha = 0.2, x = 1.0, y = 0.05, ysize = 0.1) +
        add_phylopic(agoutiPic, alpha = 0.2, x = 2.0, y = 0.05, ysize = 0.1)

    ggsave(
        filename = "../Figures/MultispeciesModeling/CommunityPredictions.png",
        plot = p,
        width = 8, height = 4
    )
}




################################################################################
############### PLOTTING CONDITIONAL OCCUPANCY FOR ALL SPECIES #################
################################################################################


# PLOT BY COMMUNITY WITH THE COMMUNITY-ONLY MODEL
# make a dataframe with column of species 1 and column of species 2 with all combinations of species
communities <- c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
community_conditional_occupancy <- data.frame(
    Species1 = rep(expand.grid(commonNames, commonNames)$Var2, each = length(communities)),
    Species2 = rep(expand.grid(commonNames, commonNames)$Var1, each = length(communities)),
    Occupancy = NA,
    SE = NA,
    lower = NA,
    upper = NA
)
community_conditional_occupancy <- community_conditional_occupancy[community_conditional_occupancy$Species1 != community_conditional_occupancy$Species2, ]
community_conditional_occupancy$Community <- rep(communities, times = nrow(community_conditional_occupancy) / length(communities))

# combine these the long way
commPresent <- community_conditional_occupancy
commAbsent <- community_conditional_occupancy
commPresent$conditionalLabel <- paste0(commPresent$Species2," present")
commAbsent$conditionalLabel <- paste0(commPresent$Species2, " absent")
commPresent$conditional <- "Species2 Present"
commAbsent$conditional <- "Species2 Absent"
community_conditional_occupancy <- rbind(commPresent, commAbsent)




# prediction data frame
dfCommunity <- dfTemplate
dfCommunity$Community <- rep(communities,
    each = 100 / length(unique(communities))
)
dfCommunity <- dfCommunity[!duplicated(dfCommunity), ] # simplify
speciesDictionary <- data.frame(
    CasualName = casualNames,
    CommonName = commonNames
)


# fill in community_conditional_occupancy
for (i in 1:nrow(community_conditional_occupancy)) {
    # establish variables
    Species1 <- speciesDictionary$CasualName[speciesDictionary$CommonName == community_conditional_occupancy$Species1[i]]
    Species2 <- speciesDictionary$CasualName[speciesDictionary$CommonName == community_conditional_occupancy$Species2[i]]
    communityInQuestion <- community_conditional_occupancy$Community[i]

    if (community_conditional_occupancy$conditional[i] == "Species2 Present") {
        # prediction where both species present
        pred <- c() # resets it so no accidental carryover
        pred <- predict(community_mod_penalty,
            type = "state",
            species = Species1,
            cond = Species2,
            newdata = dfCommunity
        )
        pred$Community <- dfCommunity$Community
        community_conditional_occupancy$Occupancy[i] <- pred$Predicted[pred$Community == communityInQuestion]
        community_conditional_occupancy$SE[i] <- pred$SE[pred$Community == communityInQuestion]
        community_conditional_occupancy$lower[i] <- pred$lower[pred$Community == communityInQuestion]
        community_conditional_occupancy$upper[i] <- pred$upper[pred$Community == communityInQuestion]
    } else if (community_conditional_occupancy$conditional[i] == "Species2 Absent") {
        # prediction where Species2 is absent
        pred <- c() # resets it so no accidental carryover
        pred <- predict(community_mod_penalty,
            type = "state",
            species = Species1,
            cond = paste0("-", Species2),
            newdata = dfCommunity
        )
        pred$Community <- dfCommunity$Community
        community_conditional_occupancy$Occupancy[i] <- pred$Predicted[pred$Community == communityInQuestion]
        community_conditional_occupancy$SE[i] <- pred$SE[pred$Community == communityInQuestion]
        community_conditional_occupancy$lower[i] <- pred$lower[pred$Community == communityInQuestion]
        community_conditional_occupancy$upper[i] <- pred$upper[pred$Community == communityInQuestion]
    } else 
    next
}

# formatting to prep for plotting
community_conditional_occupancy$Community <- gsub("Zabalo", "Zábalo", x = community_conditional_occupancy$Community)
community_conditional_occupancy$Community <- factor(community_conditional_occupancy$Community,
    levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)
community_conditional_occupancy$Species1 <- factor(community_conditional_occupancy$Species1, levels = commonNames) # so plotting doesn't alphabetize species
community_conditional_occupancy$Species2 <- factor(community_conditional_occupancy$Species2, levels = commonNames) # so plotting doesn't alphabetize species

# make labels per covariate
colors <- c(
    "Zábalo" = "darkgreen", "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen", "San Pablo" = "gold1", "Siona" = "darkgoldenrod3"
)

# make titles per potential species
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


# plot interactions for each species
for(i in 1:length(casualNames)){
    dataSub <- subset(community_conditional_occupancy, Species1 == commonNames[i])
    possibleInteractions <- unique(dataSub$Species2)

    interactionPlots <- list()
    for(j in 1:length(possibleInteractions)){
    # subset the data again
    dataSubSub <- subset(dataSub, Species2 == possibleInteractions[j])


    # facet by Species2
    dodge <- position_dodge(width = 0.3)
    p <- ggplot(dataSub, aes(x = Community, y = Occupancy, color = conditionalLabel)) +
    facet_wrap(~Species2) + 
    geom_point(aes(color = conditionalLabel), position = dodge, size = 1.5) +
    geom_errorbar(
        aes(
            ymin = lower,
            ymax = upper,
            color = conditionalLabel
        ),
        position = dodge, width = 0.15, linewidth = .5
    ) +
    labs(x = "Community", 
    y = paste0("Conditional ", tolower(commonNames[i]), " occupancy probability")) +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(
        text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)
    ) 
    interactionPlots[[j]] <- p
}

# make the plot
n <- length(interactionPlots)
nCol <- floor(sqrt(n))
p <- do.call("grid.arrange", c(interactionPlots, ncol = nCol))
annotate_figure(p, top = text_grob(paste0(commonNames[i], " Interactions (Penalized community-only model)"),
    face = "bold", size = 14
))
ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[i], "_InteractionsByCommunity.png"),
    width = 10, height = 8
)
}








# TESTING
# use the predict function without providing new data to predict by site, then extrapolate to by community
test <- as.data.frame(predict(natural_mod_penalty,
    type = "state",
    species = "agouti",
    cond = "paca"
))
test$Site <- rownames(combinedDetection)

test$Community <- NA
for(i in 1:nrow(test)){
    # if a site starts with XYZ, then it's in the XYZ community
    if (grepl("SNA", test$Site[i])){
        test$Community[i] <- "Siona"
    } else if(grepl("ZAB", test$Site[i])){
        test$Community[i] <- "Zábalo"
    } else if(grepl("SGE", test$Site[i])){
        test$Community[i] <- "Sinangoe"
    } else if(grepl("SKP", test$Site[i])){
        test$Community[i] <- "Remolino"
    }

    # now assign San Pablo: SKP31-37, SNA3
    if (test$Site[i] == "SNA3" | test$Site[i] == "SKP31" | test$Site[i] == "SKP32" |
        test$Site[i] == "SKP33" | test$Site[i] == "SKP34" | test$Site[i] == "SKP35" |
        test$Site[i] == "SKP36" | test$Site[i] == "SKP37") {
        test$Community[i] <- "San Pablo"
    }
}

# group by community and average the occupancy and calculate 95% CI
test2 <- test %>%
    group_by(Community) %>%
    summarize(
        SD = sd(Predicted),
        Predicted = mean(Predicted), 
        count = n()) %>%
    mutate(
        SE = SD / sqrt(count),
        lower = Predicted - qt(0.975, df = count - 1) * SD / sqrt(count),
        upper = Predicted + qt(0.975, df = count - 1) * SD / sqrt(count)
    )
test2


# OPTION TWO: use the predict function with new data to predict by community






