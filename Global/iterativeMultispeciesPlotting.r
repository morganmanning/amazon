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
require(cowplot) 
require(grid) 

rm(list = ls())

# define the species pairs to iterate through
species_pairs <- c(
    "peccary_agouti_MSM.RData",
    "paca_agouti_MSM.RData",
    "ocelot_agouti_MSM.RData"
)

# store all interaction plots for the combined figure
all_interaction_plots <- list()
plot_counter <- 1

# magnitude of difference
difference <- read.csv("hedges_g.csv")
difference$Community <- factor(difference$Community,
    levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
)

# colors outside the loop since they're the same
colors <- c(
    "Zábalo" = "darkgreen",
    "Remolino" = "forestgreen",
    "Sinangoe" = "yellowgreen",
    "San Pablo" = "gold1",
    "Siona" = "darkgoldenrod3"
)

# make axis titles per potential species
peccary <- ~ atop(paste("Collared peccary"), paste("(", italic("Pecari tajacu"), ")"))
paca <- ~ atop(paste("Lowland paca"), paste("(", italic("Cuniculus paca"), ")"))
agouti <- ~ atop(paste("Black agouti"), paste("(", italic("Dasyprocta fuliginosa"), ")"))
ocelot <- ~ atop(paste("Ocelot"), paste("(", italic("Leopardus pardalis"), ")"))

peccaryPic <- get_uuid(name = "Pecari tajacu", n = 1)
pacaPic <- get_uuid(name = "Cuniculus paca", n = 1)
agoutiPic <- get_uuid(name = "Dasyprocta", n = 1)
ocelotPic <- get_uuid(name = "Leopardus pardalis", n = 1)

# iterate through each species pair
for (pair_file in species_pairs) {
    # load data (the objects inside have the same names)
    # use a temporary environment to avoid conflicts
    temp_env <- new.env()
    load(paste0("R Objects/", pair_file), envir = temp_env)

    # get the objects from the temporary environment
    casualNames <- temp_env$casualNames
    commonNames <- temp_env$commonNames
    latinNames <- temp_env$species
    null_multispecies_model <- temp_env$null_multispecies_model
    community_mod_penalty <- temp_env$community_mod_penalty

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

    ggsave(paste0("../Figures/MultispeciesModeling/", paste0(casualNames, collapse = ""), "NullModelPredictions.png"),
        width = 5, height = 5
    )

    ################################################################################
    ############################### PLOT EFFECT SIZES ##############################
    ################################################################################

    sum_out <- summary(community_mod_penalty)

    # helper function to convert coefficient matrix into a tidy data frame
    tidy_coefs <- function(mat, type) {
        as.data.frame(mat) |>
            tibble::rownames_to_column("term") |>
            dplyr::mutate(
                param_type = type,
                estimate = Estimate,
                lower = Estimate - 1.96 * SE,
                upper = Estimate + 1.96 * SE,
                Species = gsub("^\\[|\\].*", "", term),
                param = gsub(".*\\] ", "", term),
                significant = ifelse(`P(>|z|)` < 0.05, "Significant", "Not Significant")
            ) |>
            dplyr::select(Species, param, estimate, lower, upper, significant, param_type)
    }

    # tidy up both occupancy and detection
    occ_df <- tidy_coefs(sum_out$state, "Occupancy")
    det_df <- tidy_coefs(sum_out$det, "Detection")

    # combine them and make zabalo the intercept
    longDF <- dplyr::bind_rows(occ_df, det_df)
    longDF$param <- gsub("^\\(Intercept\\)$", "CommunityZabalo", longDF$param)

    # just occupancy
    longDFOccupancy <- longDF[longDF$param_type == "Occupancy", ]

    # now plot it
    ggplot(longDFOccupancy, aes(x = estimate, y = Species, color = significant)) +
        geom_point(size = 2) +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        facet_wrap(~param, labeller = label_wrap_gen(width = 20)) +
        labs(x = "Effect size (95% CI)", y = NULL) +
        scale_color_manual(values = c("Significant" = "darkred", "Not Significant" = "gray40")) +
        theme_minimal() +
        theme(
            text = element_text(family = "Times", colour = "black"),
            legend.title = element_blank(),
            strip.text = element_text(size = 10),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 12),
            panel.border = element_rect(color = "black", size = 0.5, fill = "transparent")
        )

    # save it
    ggsave(paste0("../Figures/MultispeciesModeling/", paste0(casualNames, collapse = ""), "CommunityPenalizedEffectSizes.png"),
        width = 7, height = 5
    )

    ##### TABLE
    # prep the table
    table_df <- longDF %>%
        filter(param_type == "Occupancy") %>%
        mutate(
            Covariate = param,
            `95% CI` = ifelse(
                !is.na(lower) & !is.na(upper),
                sprintf("[%.2f, %.2f]", lower, upper),
                "-"
            ),
            Estimate_fmt = ifelse(
                !is.na(lower) & !is.na(upper) & (lower > 0 | upper < 0),
                sprintf("<b>%.3f</b>", estimate),
                sprintf("%.3f", estimate)
            )
        ) %>%
        arrange(Species, Covariate) %>%
        group_by(Species) %>%
        mutate(
            Species_display = ifelse(row_number() == 1, Species, "")
        ) %>%
        ungroup() %>%
        select(
            Species = Species_display,
            Covariate,
            Estimate = Estimate_fmt,
            `95% CI`
        )

    kbl(table_df,
        col.names = c("Species", "Covariate", "Estimate", "95% CI"),
        escape = FALSE
    ) %>%
        kable_classic(full_width = TRUE, html_font = "Times New Roman") %>%
        column_spec(1, bold = TRUE, italic = TRUE) %>%
        save_kable(paste0("../Figures/MultispeciesModeling/", paste0(casualNames, collapse = ""), "CommunityPenalizedEffectSizesTable.png"),
            zoom = 1.5
        )

    ################################################################################
    #################### CONDITIONAL OCCUPANCY PREDICTIONS #########################
    ################################################################################

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

    # natural_conditional_occupancy <- data.frame(
    #     Species1 = expand.grid(casualNames, casualNames)$Var2,
    #     Species2 = expand.grid(casualNames, casualNames)$Var1,
    #     Present1_Present2 = NA,
    #     PresentSE = NA,
    #     PresentLower = NA,
    #     PresentUpper = NA,
    #     Present1_Absent2 = NA,
    #     AbsentSE = NA,
    #     AbsentLower = NA,
    #     AbsentUpper = NA
    # )
    # natural_conditional_occupancy <- natural_conditional_occupancy[natural_conditional_occupancy$Species1 != natural_conditional_occupancy$Species2, ]

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

    # fill in conditional occupancy predictions for null model
    for (i in 1:nrow(conditional_occupancy)) {
        # prediction where Species2 is present
        pred <- predict(null_multispecies_model,
            type = "state",
            species = conditional_occupancy$Species1[i],
            cond = conditional_occupancy$Species2[i]
        )[1, ] # pull just the first row since rows are duplicates
        conditional_occupancy$Present1_Present2[i] <- pred$Predicted
        conditional_occupancy$PresentSE[i] <- pred$SE
        conditional_occupancy$PresentLower[i] <- pred$lower
        conditional_occupancy$PresentUpper[i] <- pred$upper

        # prediction where Species2 is absent
        pred <- predict(null_multispecies_model,
            type = "state",
            species = conditional_occupancy$Species1[i],
            cond = paste0("-", conditional_occupancy$Species2[i])
        )[1, ] # pull just the first row since rows are duplicates
        conditional_occupancy$Present1_Absent2[i] <- pred$Predicted
        conditional_occupancy$AbsentSE[i] <- pred$SE
        conditional_occupancy$AbsentLower[i] <- pred$lower
        conditional_occupancy$AbsentUpper[i] <- pred$upper
    }


    ################################################################################
    ################ CONDITIONAL OCCUPANCY BY COMMUNITY ############################
    ################################################################################

    # create new data frame with community as covariate
    dfCommunity <- data.frame(
        Community = c("Zabalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    )

    community_conditional_occupancy <- data.frame(
        Species1 = rep(expand.grid(commonNames, commonNames)$Var2, each = length(dfCommunity$Community)),
        Species2 = rep(expand.grid(commonNames, commonNames)$Var1, each = length(dfCommunity$Community)),
        Occupancy = NA,
        SE = NA,
        lower = NA,
        upper = NA
    )
    community_conditional_occupancy <- community_conditional_occupancy[community_conditional_occupancy$Species1 != community_conditional_occupancy$Species2, ]
    community_conditional_occupancy$Community <- rep(dfCommunity$Community, times = nrow(community_conditional_occupancy) / length(dfCommunity$Community))

    # combine these the long way
    commPresent <- community_conditional_occupancy
    commAbsent <- community_conditional_occupancy
    commPresent$conditionalLabel <- paste0(commPresent$Species2, " present")
    commAbsent$conditionalLabel <- paste0(commPresent$Species2, " absent")
    commPresent$conditional <- "Species2 Present"
    commAbsent$conditional <- "Species2 Absent"
    community_conditional_occupancy <- rbind(commPresent, commAbsent)

    # loop through each community and species combination
    for (i in 1:nrow(community_conditional_occupancy)) {
        communityInQuestion <- community_conditional_occupancy$Community[i]
        Species1 <- community_conditional_occupancy$Species1[i]
        Species2 <- community_conditional_occupancy$Species2[i]

        # match common names to casual names
        Species1_casual <- casualNames[commonNames == Species1]
        Species2_casual <- casualNames[commonNames == Species2]

        if (community_conditional_occupancy$conditional[i] == "Species2 Present") {
            # prediction where Species2 is present
            pred <- predict(community_mod_penalty,
                type = "state",
                species = Species1_casual,
                cond = Species2_casual,
                newdata = dfCommunity
            )
            pred$Community <- dfCommunity$Community
            community_conditional_occupancy$Occupancy[i] <- pred$Predicted[pred$Community == communityInQuestion]
            community_conditional_occupancy$SE[i] <- pred$SE[pred$Community == communityInQuestion]
            community_conditional_occupancy$lower[i] <- pred$lower[pred$Community == communityInQuestion]
            community_conditional_occupancy$upper[i] <- pred$upper[pred$Community == communityInQuestion]
        } else if (community_conditional_occupancy$conditional[i] == "Species2 Absent") {
            # prediction where Species2 is absent
            pred <- predict(community_mod_penalty,
                type = "state",
                species = Species1_casual,
                cond = paste0("-", Species2_casual),
                newdata = dfCommunity
            )
            pred$Community <- dfCommunity$Community
            community_conditional_occupancy$Occupancy[i] <- pred$Predicted[pred$Community == communityInQuestion]
            community_conditional_occupancy$SE[i] <- pred$SE[pred$Community == communityInQuestion]
            community_conditional_occupancy$lower[i] <- pred$lower[pred$Community == communityInQuestion]
            community_conditional_occupancy$upper[i] <- pred$upper[pred$Community == communityInQuestion]
        }
    }

    # add a more descriptive label for plotting
    community_conditional_occupancy$conditionalLabel <- ifelse(
        community_conditional_occupancy$conditional == "Species2 Present",
        paste0(community_conditional_occupancy$Species2, " present"),
        paste0(community_conditional_occupancy$Species2, " absent")
    )

    # formatting to prep for plotting
    community_conditional_occupancy$Community <- gsub("Zabalo", "Zábalo", x = community_conditional_occupancy$Community)
    community_conditional_occupancy$Community <- factor(community_conditional_occupancy$Community,
        levels = c("Zábalo", "Remolino", "Sinangoe", "San Pablo", "Siona")
    )
    community_conditional_occupancy$Species1 <- factor(community_conditional_occupancy$Species1, levels = commonNames)
    community_conditional_occupancy$Species2 <- factor(community_conditional_occupancy$Species2, levels = commonNames)

    # save it to facilitate plotting
    #save(community_conditional_occupancy,
        #file = paste0("R Objects/", gsub(" ", "", casualNames[1]), gsub(" ", "", casualNames[2]), "MSM_byCommunity.RData")
    #)


    # plot interactions for each species
    for (i in 1:length(casualNames)) {
        dataSub <- subset(community_conditional_occupancy, Species1 == commonNames[i])
        possibleInteractions <- unique(dataSub$Species2)
        possibleInteractionLatinName <- latinNames[commonNames == possibleInteractions]

        for (j in 1:length(possibleInteractions)) {
            # subset the data again
            dataSubSub <- subset(dataSub, Species2 == possibleInteractions[j])

            dodge <- position_dodge(width = 0.3)

            # create the interaction plot
            p <- ggplot(dataSubSub, aes(x = Community, y = Occupancy, color = Community, linetype = conditionalLabel)) +
                geom_point(aes(color = Community), position = dodge, size = 1.5) +
                geom_errorbar(
                    aes(
                        ymin = lower,
                        ymax = upper,
                        color = Community,
                        linetype = conditionalLabel
                    ),
                    position = dodge,
                    width = 0.15,
                    linewidth = .5
                ) +
                labs(
                    x = "Territory",
                    y = "Conditional occupancy probability",
                    title = paste0(commonNames[i], " | ", possibleInteractions[j])
                ) +
                ylim(c(0, 1)) +
                theme_classic() +
                scale_color_manual(values = colors) +
                scale_linetype_manual(values = c("dashed", "solid")) +
                theme(
                    text = element_text(family = "Times", colour = "black"),
                    axis.text = element_text(colour = "black"),
                    legend.title = element_blank(),
                    legend.position = "top",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2),
                    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
                ) +
                guides(color = "none")

            # phylopic additions
            p <- p +
                add_phylopic(uuid = ifelse(commonNames[i] == "Black agouti",
                    get_uuid(name = "Dasyprocta", n = 1),
                    get_uuid(name = latinNames[i], n = 1)
                ), alpha = 0.2, x = 1.2, y = 0.05, ysize = 0.1) +
                add_phylopic(uuid = ifelse(possibleInteractions[j] == "Black agouti",
                    get_uuid(name = "Dasyprocta", n = 1),
                    get_uuid(name = possibleInteractionLatinName[j], n = 1)
                ), alpha = 0.2, x = 2.1, y = 0.05, ysize = 0.1)
            # add a vertical line to symbolize conditionality (idk if I like this)
            p <- p +
                geom_segment(aes(x = 1.65, xend = 1.65, y = 0.001, yend = 0.1),
                    color = "black",
                    size = 0.2
                )

            # store the plot for the combined figure
            all_interaction_plots[[plot_counter]] <- p
            plot_counter <- plot_counter + 1

            # also save individual plots as before (only for two species pairs)
            if (length(casualNames) == 2) {
                otherSpecies <- casualNames[casualNames != casualNames[i]]

                # create individual plot with grid.arrange
                individual_plot <- arrangeGrob(p)

                ggsave(paste0("../Figures/MultispeciesModeling/", casualNames[i], "_", otherSpecies, "_InteractionsByCommunity.png"),
                    plot = individual_plot,
                    width = 10, height = 6
                )
            }
        }
    }


    #################### PLOT DIFFERENCE IN MEANS WITH SE FOR EACH SPECIES PAIR
    differenceSpecies <- unique(difference$Species1)
    communities <- unique(difference$Community)

    # prevent issues
    difference$Species1 <- factor(difference$Species1, levels = differenceSpecies)
    difference$Species2 <- factor(difference$Species2, levels = differenceSpecies)
    differenceSpecies <- factor(differenceSpecies, levels = differenceSpecies)

    # make a plot for each species pair showing difference in means at each community
    for (i in 1:length(differenceSpecies)) {
        dataSub <- subset(difference, Species1 == differenceSpecies[i])
        possibleInteractions <- unique(dataSub$Species2)

        for (j in 1:length(possibleInteractions)) {
            dataSubSub <- subset(dataSub, Species2 == possibleInteractions[j])

            dodge <- position_dodge(width = 0.3)

            p <- ggplot(dataSubSub, aes(x = Community, y = meanDifference, color = Community)) +
                geom_point(aes(color = Community), position = dodge, size = 1.5) +
                geom_errorbar(
                    aes(
                        ymin = meanDifference - meanDifferenceSE,
                        ymax = meanDifference + meanDifferenceSE,
                        color = Community
                    ),
                    position = dodge,
                    width = 0.15,
                    linewidth = .5
                ) +
                labs(
                    x = "Territory",
                    y = paste0(
                        "Difference in mean conditional ", tolower(differenceSpecies[i]),
                        " occupancy probability"
                    )
                ) +
                ylim(c(-0.1, 1)) +
                scale_color_manual(values = colors) +
                theme_classic() +
                theme(
                    text = element_text(family = "Times", colour = "black"),
                    axis.text = element_text(colour = "black"),
                    legend.title = element_blank(),
                    legend.position = "top",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major.y = element_line(color = "#cecece", linewidth = 0.2)
                )

            p_with_title <- annotate_figure(p,
                top = text_grob(
                    paste0(
                        "Difference in mean conditional ", tolower(differenceSpecies[i]),
                        " occupancy when ", tolower(possibleInteractions[j]), " is present vs. absent \n",
                        "(Penalized territory-only model)"
                    ),
                    face = "bold",
                    size = 14
                )
            )

            figurePath <- paste0(
                "../Figures/MultispeciesModeling/",
                gsub(" ", "", tolower(differenceSpecies[i])), "_",
                gsub(" ", "", tolower(possibleInteractions[j])), "_DifferencesByCommunity.png"
            )

            ggsave(figurePath, plot = p_with_title, width = 10, height = 8)
        }
    }
} # main loop through species pairs

################################################################################
################ CREATE COMBINED FACETED PLOT #################################
################################################################################

# # ensure all plots have the same y-axis label for consistency
# for (i in 1:length(all_interaction_plots)) {
#     all_interaction_plots[[i]] <- all_interaction_plots[[i]] +
#         labs(y = "Conditional occupancy probability")
# }

# create the combined 3x2 faceted plot with shared y-axis
combined_plot <- plot_grid(
    plotlist = all_interaction_plots,
    ncol = 2,
    nrow = 3,
    align = "hv", # align both horizontally and vertically
    axis = "lr", # share left-right axes
    rel_heights = c(1, 1, 1),
    rel_widths = c(1, 1),
    labels = c("A", "B", "C", "D", "E", "F"),
    label_size = 14,
    label_fontfamily = "Times"
)

# add overall title and axis labels
combined_plot_labelled <- ggdraw() +
    draw_plot(combined_plot, 0.025, 0, 0.97, 1) +
    # draw_label(
    #     "All Species Interactions by Territory",
    #     x = 0.5, y = 0.98,
    #     size = 16,
    #     fontfamily = "Times",
    #     fontface = "bold"
    # ) +
    draw_label(
        "Conditional occupancy probability",
        x = 0.01, y = 0.5,
        size = 16,
        fontfamily = "Times",
        angle = 90
    )
combined_plot_labelled
# save it
ggsave(
    "../Figures/MultispeciesModeling/allSpeciesInteractionsByCommunity.png",
    plot = combined_plot_labelled,
    width = 8,
    height = 10,
    dpi = 300
)
