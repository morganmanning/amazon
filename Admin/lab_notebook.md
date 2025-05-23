## Methods Briefly

## Single species models
### Detection matrices
- For each species, the null model was run grouping detection into occasions spanning 1 day (i.e. ungrouped) to 20 days. 
    - The optimal number of days per occasion was chosen based on the null model with the lowest occupancy estimate standard error. 
    - The grouping number can be found per species within occupancyEstimatesPlots.R in the object 'nDaysGroupedPerSpecies' (also included in occupancyDetectionEstimates.png)
    - The optimal grouping number was used to consolidate a detection history for each species and used moving forward for occupancy modeling

### Detection covariates
- Community and EffortScaled were considered as covariates of detection for all species
    - Ran models with both together, each individually, then the null (~1) with occupancy having no covariates (~1)
    - The best model (per AIC) for detection was then used for consideration of occupancy covariates

### Occupancy covariates
- Using the best detection model per species, the following covariates were considered for occupancy: the community where the camera trap is associated with, the percent natural area within 25 km (https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d; ESRI Sentinel-2 10m time series; Coordinate system: WGS84; UTM; ESPG:3857; water, trees, and flodded vegetations were considered "natural area"), the amount of rainfall at the camera trap location (scaled; https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary), the distance to a water source (derived from https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d), the temperature at the camera trap location (scaled; https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary), and the distance to the nearest community.
    - All combinations of these covariates (including the null) were considered for each species.
    - Models with convergence issues were removed from consideration
    - All best models (within 2 AIC of the best model) were averaged 
- Null models were used to estimate naive detection and occupancy for each species (there's a table of and plot)
- The best models averaged were used to back-predict occupancy assuming the average of each covariate and then that occupancy was averaged to get an estimate of occupancy considering existing conditions/covariates
- The best models averaged were used to predict species occupancy across each covariate at each community


## Multispecies models
- Species interactions looked at:
    + predator/prey: ocelot/agouti
    + competitors: agouti/paca
    + two hunted species but peccary is more desirable: agouti/peccary
- The species type interactions were chosen to answer ecological/biological questions about how animal interactions are impacted across communities, but the specific species were chosen because they had the highest number of detections within the species type
- Detection histories were condensed to where each trapping occasion equaled two days
    + this was picked because it was the most stable grouping across the study species and detection histories needed to have the same dimensions (i.e., I had to pick one grouping for all species)
- Ran multispecies models for each pairing using the unmarked package in R
- Due to low detection, the inclusion of environmental covariates was ruled out. When an environmental covariate was included, the model broke.
    + Instead, only the community was included as a predictor of both occupancy and detection. 
- Utilized penalized likelihood with the community model to improve parameter estimates (https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210X.12368).
    + Identified the optimal value of the bayes-inspired penalty term for the models that support penalized likelihood. 
    + For each potential value of the penalty term (0.5 and 1), K-fold cross validation is performed.
    + Log-likelihoods for the test data in each fold are calculated and summed. 
    + The penalty term that maximizes the sum of the fold log-likelihoods is selected as the optimal value. 
    + Finally, the community model is re-fit with the full dataset using the selected penalty term. 
- Using the penalized model, occupancy of each species was predicted across all communities in two scenarios: where the paired species was present and where the paired species was absent
- The magnitude of difference in occupancy predictions was quantified for each species pairing at each community using the difference in means with the standard error calculated accounting for the unequal variances among communities (https://math.stackexchange.com/questions/1257855/two-different-formulas-for-standard-error-of-difference-between-two-means)


### Multidimensional scaling
- Because the data did not allow for the inclusion of covariates in the multispecies models, we opted to use multivariate analysis, specifically non-metric multidimensional scaling, to examine the relationship among communities. 
- The NMDS results will help us better understand the differences in wildlife interactions among communities, as a result of community characteristics.
- The covariates that were included in the NMDS include geographical, ecological, and sociological forces:
    + The distance from the community to the nearest community
    + The area of the indigenous territory
    + Rainfall (Averaged within the community/territory between July 2020 and December 2020 (the months that the cameras were out); https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary)
    + Humidity (Averaged within the community/territory between July 2020 and December 2020 (the months that the cameras were out); https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary)
    + Air temperature (Averaged within the community/territory between July 2020 and December 2020 (the months that the cameras were out); https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary)
    + Root moisture (Averaged within the community/territory between July 2020 and December 2020 (the months that the cameras were out); https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary)
    + Average percent natural area within 25 km of camera traps (https://www.arcgis.com/home/item.html?id=cfcb7609de5f478eb7666240902d4d3d; ESRI Sentinel-2 10m time series; Coordinate system: WGS84; UTM; ESPG:3857; water, trees, and flodded vegetations were considered "natural area")
    + The total number of individual detections
    + The number of species detected
    + Shannon diversity index
    + Simpson diversity index
    + Population size (don't have data for San Pablo and Remolino)
    + Lat and lon
    + Mean elevation within territory (https://www.sciencebase.gov/catalog/item/5920dd83e4b0ac16dbdf3a4d)
    + Potentials: Average number of days spent hunting per person per month (divided by wet and dry season; verify this), average number of days spent fishing per person per month (divided by wet and dry season; verify this), percent of the population who hunt, percent of the population who fish, lat and lon, standard deviation of root moisture/air temp/humidity/rainfall
- The a distance matrix was calculated pairwise for all covariates, then we used multidimensional scaling
- used "vegan" package in R to do NMDS



## Other ideas:
- sensitivity analysis with percent natural area looking at different buffer distances
- plot null occupancies versus average occupancy considering covariates on the same plot?
- move everything to Bayesian instead of unmarked
- find better Ecuador shapefile for plotting sites
    - plot territories on site map?
-