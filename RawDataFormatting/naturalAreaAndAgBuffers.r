# set up
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# wd
setwd("Documents/amazon/Global/Data/Community-level Covariates") 

######## SOURCE
# data for this came from the sentinel-2 10m land use/land cover time series dataset
# natural area = water, trees, flooded vegetation, snow/ice
# ag = crops, rangeland

################################################################################
########################## PROCESS AGRICULTURE #################################
################################################################################
# all agriculture CSVs
ag_files <- list.files(pattern = "^\\d+KM_ag\\.csv$")

# process each file
ag_list <- list()
for (file in ag_files) {
    # extract distance from filename
    distance <- str_extract(file, "^\\d+KM")

    # read CSV
    df <- read_csv(file, show_col_types = FALSE)

    # calculate proportion of area
    sum_col <- names(df)[str_detect(names(df), "_sum")]
    count_col <- names(df)[str_detect(names(df), "_count")]

    df <- df %>%
        mutate(proportion = !!sym(sum_col) / !!sym(count_col)) %>%
        select(!!sym(names(df)[1]), proportion) # keep ID + proportion

    # rename proportion column to distance of buffer
    colnames(df)[2] <- distance

    # store it
    ag_list[[distance]] <- df
}

# merge 
ag_table <- reduce(ag_list, full_join, by = names(ag_list[[1]])[1])

# save
write_csv(ag_table, "agricultureBuffers.csv")


################################################################################
########################## PROCESS NATURAL AREA ################################
################################################################################ 
nat_files <- list.files(pattern = "^\\d+KM_natArea\\.csv$")

nat_list <- list()
for (file in nat_files) {
    # extract distance from filename
    distance <- str_extract(file, "^\\d+KM")

    # read CSV
    df <- read_csv(file, show_col_types = FALSE)

    # calculate proportion of area
    sum_col <- names(df)[str_detect(names(df), "_sum")]
    count_col <- names(df)[str_detect(names(df), "_count")]

    # calculate proportion of area
    df <- df %>%
        mutate(proportion = !!sym(sum_col) / !!sym(count_col)) %>%
        select(!!sym(names(df)[1]), proportion)

    # rename proportion column to distance of buffer
    colnames(df)[2] <- distance

    # store it
    nat_list[[distance]] <- df
}

nat_table <- reduce(nat_list, full_join, by = names(nat_list[[1]])[1])
write_csv(nat_table, "naturalAreaBuffers.csv")



################################################################################
########################## SENSITIVITY ANALYSIS ################################
################################################################################

# read in data
ag_table <- read_csv("agricultureBuffers.csv")
nat_table <- read_csv("naturalAreaBuffers.csv")

# check correlation between buffers
ag_cor <- cor(ag_table[, -1], use = "complete.obs")
nat_cor <- cor(nat_table[, -1], use = "complete.obs")

print("Agriculture buffer correlations:")
print(round(ag_cor, 3))
print("\nNatural area buffer correlations:")
print(round(nat_cor, 3))

# coefficient of variation across buffers (site-level)
calc_cv <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100

ag_cv <- apply(ag_table[, -1], 1, calc_cv)
nat_cv <- apply(nat_table[, -1], 1, calc_cv)

print(paste("\nMean CV across buffers - Agriculture:", round(mean(ag_cv), 2), "%"))
print(paste("Mean CV across buffers - Natural area:", round(mean(nat_cv), 2), "%"))

# variation between sites at each buffer
ag_var <- apply(ag_table[, -1], 2, function(x) sd(x, na.rm = TRUE))
nat_var <- apply(nat_table[, -1], 2, function(x) sd(x, na.rm = TRUE))

print("\nStandard deviation between sites:")
print("Agriculture:")
print(round(ag_var, 3))
print("Natural area:")
print(round(nat_var, 3))

##### PLOT IT
library(ggplot2)
library(tidyr)

# reshape for plotting
ag_long <- ag_table %>%
    pivot_longer(-1, names_to = "buffer", values_to = "proportion") %>%
    mutate(buffer_km = as.numeric(str_extract(buffer, "\\d+")))

nat_long <- nat_table %>%
    pivot_longer(-1, names_to = "buffer", values_to = "proportion") %>%
    mutate(buffer_km = as.numeric(str_extract(buffer, "\\d+")))

# plot
p1 <- ggplot(ag_long, aes(x = buffer_km, y = proportion, group = 1, color = factor(1))) +
    geom_line() +
    geom_point() +
    facet_wrap(~ get(names(ag_table)[1]), ncol = 1) +
    labs(title = "Agriculture by Buffer Size", x = "Buffer (km)", y = "Proportion") +
    theme_minimal() +
    theme(legend.position = "none")

p2 <- ggplot(nat_long, aes(x = buffer_km, y = proportion, group = 1, color = factor(1))) +
    geom_line() +
    geom_point() +
    facet_wrap(~ get(names(nat_table)[1]), ncol = 1) +
    labs(title = "Natural Area by Buffer Size", x = "Buffer (km)", y = "Proportion") +
    theme_minimal() +
    theme(legend.position = "none")

print(p1)
print(p2)

# Correlations (>0.9): If adjacent buffers are highly correlated, they're providing similar information. You could use a smaller buffer to save computation time.
# Correlations (<0.7): If larger buffers show weak correlations with smaller ones, buffer choice matters significantly and you should justify your selection based on ecological theory.
# CV values: High CV (>20%) means proportions change substantially across buffer sizes for individual sites, suggesting buffer choice is important. Low CV means results are robust to buffer selection.
# Between-site variation: The buffer that maximizes variation between sites (highest SD) may be most useful for detecting differences among communities.
# Visual patterns: Look for sites where proportion plateaus at a certain buffer size - that suggests the relevant spatial scale for that community.