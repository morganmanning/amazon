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
- Using the best detection model per species, the following covariates were considered for occupancy: the community where the camera trap is associated with, the percent natural area within 25 km, the amount of rainfall at the camera trap location (scaled), the distance to a water source, the temperature at the camera trap location (scaled), and the distance to the nearest community.
    - All combinations of these covariates (including the null) were considered for each species.
    - Models with convergence issues were removed from consideration
    - All best models (within 2 AIC of the best model) were averaged 
- Null models were used to estimate naive detection and occupancy for each species (there's a table of and plot)
- The best models averaged were used to back-predict occupancy assuming the average of each covariate and then that occupancy was averaged to get an estimate of occupancy considering existing conditions/covariates
- The best models averaged were used to predict species occupancy across each covariate at each community

## Multispecies models


## Other ideas:
- sensitivity analysis with percent natural area looking at different buffer distances
- plot null occupancies versus average occupancy considering covariates on the same plot?
- move everything to Bayesian instead of unmarked