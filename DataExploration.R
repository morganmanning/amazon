setwd("~/Documents/amazon/Source Material")

require(dplyr)

###################################################
############### COLLARED PECCARY ##################
###################################################

# load data
species <- read.csv("CollaredPeccary.csv")
species <- data.frame(species[,-1], row.names=species[,1]) #rownames = stations

# occupancy 
library(unmarked)

# data as matrix
y <- as.matrix(species)
y <- y[ order(as.numeric(row.names(y))), ] #order matters

# occupancy covariates
siteCovariate = read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$Community <- siteCovariate$Community/1000
siteCovariate$River <- siteCovariate$River/1000

# unmarked df
ufo = unmarkedFrameOccu(y, 
                        siteCovs = siteCovariate,
                        obsCovs = NULL)

plot(ufo)

# models (non significant)
Null = occu( ~1 ~1, ufo) # detection and occupancy are constant/unvarying
Community = occu( ~1 ~Community, ufo) # detection is constant, but occupancy varies with the distance to the community
      # backTransform(Community, type = "det") # p AKA detection
      # backTransform(Community, type = "state") # psi AKA occupancy 
      # can't do "state" because covariates are present -> would need to assign estimate for covariate(s)
River = occu( ~1 ~River, ufo) # detection is constant, but occupancy varies with the distance to the river 
Habitat <- occu(~1 ~Habitat, ufo)

CommunityRiver = occu( ~1 ~River + Community, ufo)
CommunityHabitat = occu( ~1 ~Community + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHunted = occu( ~1 ~Community + Hunting, ufo) # occupancy varies with dist to community and hunting

CommunityRiver2 = occu( ~Habitat ~River + Community, ufo)
CommunityHabitat2 = occu( ~Habitat ~Community + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHunted2 = occu( ~Habitat ~Community + Hunting, ufo) # occupancy varies with dist to community and hunting

CommunityRiverEffort = occu( ~Effort ~River + Community, ufo)
CommunityHabitatEffort = occu( ~Effort ~Community + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHuntedEffort = occu( ~Effort ~Community + Hunting, ufo) # occupancy varies with dist to community and hunting

CommunityRiverTrail = occu( ~Trail.Distance ~River + Community, ufo)
CommunityHabitatTrail = occu( ~Trail.Distance ~Community + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHuntedTrail = occu( ~Trail.Distance ~Community + Hunting, ufo) # occupancy varies with dist to community and hunting

CommunityRiver3 = occu( ~Habitat + Effort ~River + Community, ufo)
CommunityHabitat3 = occu( ~Habitat + Effort ~Community + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHunted3 = occu( ~Habitat + Effort ~Community + Hunting, ufo) # occupancy varies with dist to community and hunting

Biology <- occu(~Habitat + Effort + Trail.Distance ~ Community + Hunting + River + Trail.Distance + Habitat, ufo)

GlobalOccupancy = occu(~1 ~Community + Hunting + River + Trail.Distance + Habitat + Effort, ufo)
GlobalDetection = occu(~Community + Hunting + River + Trail.Distance + Habitat + Effort ~1, ufo)
Global <- occu(~Community + Hunting + River + Trail.Distance + Habitat + Effort 
               ~Community + Hunting + River + Trail.Distance + Habitat + Effort, 
               ufo)
#MuMIn::dredge(Global, rank = 'AIC')
# Global[1]
#str(siteCovariate)
# AIC values
BestModel = fitList(Null, Community, River, Habitat,
                    CommunityRiver, CommunityHabitat, CommunityHunted,
                    CommunityRiver2, CommunityHabitat2, CommunityHunted2,
                    CommunityRiverEffort, CommunityHabitatEffort, CommunityHuntedEffort,
                    CommunityRiverTrail, CommunityHabitatTrail, CommunityHuntedTrail,
                    CommunityRiver3, CommunityHabitat3, CommunityHunted3,
                    GlobalOccupancy, GlobalDetection, Global,
                    Biology)
modSel(BestModel)






# predict
newdata = data.frame(0:10)
colnames(newdata)[1] = "Community"
predicted = predict(Community, type="state", newdata=newdata, appendData=TRUE) # state = occupancy
par(pty="s")

# plot
library(ggplot2)
ggplot(NULL, aes(x=Community, y=Predicted)) + 
  geom_line(data=predicted, linetype="solid", size=1) +
  geom_ribbon(data=predicted, aes(ymin=lower, ymax=upper), fill="black", alpha=0.15) +
  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=1, expand = c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0:1), expand = c(0,0)) +
  theme_minimal() +
  coord_fixed(10) +
  labs(title = predicted$English, subtitle = predicted$Scientific, 
       x = "Distance from Community (km)", 
       y = "Collared Peccary Occurrence Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


levelplot(Predicted ~ ggPred$x + ggPred$y,
          data = occuPred,
          col.regions = rev(terrain.colors(100)),
          at = seq(0,1,length.out=101))


## Spatial occupancy model

stations <- read.csv("Stations2018.csv")
coords <- stations %>% 
  group_by(Station) %>%
  slice(1) %>%
  arrange(Station) %>%
  select(Station, y, x)

coords <- coords[,2:3]


library(spOccupancy)

speciesList <- list(y = species,
                    occ.covs = siteCovariate,
                    det.covs = siteCovariate,
                    coords = coords)

# # Fit a spatial, single-species occupancy model
# out.sp <- spPGOcc(occ.formula = ~ CR, 
#                   det.formula = ~ 1, 	          
#                   data = speciesList, 
#                   n.batch = 400, 
#                   batch.length = 25,
#                   n.thin = 5, 
#                   n.burn = 5000, 
#                   n.chains = 3,
#                   n.report = 100)
# summary(out.sp)




# Format with abbreviated specification of inits for alpha and beta.
inits <- list(alpha = 0, 
              beta = 0, 
              z = apply(speciesList$y, 1, max, na.rm = TRUE))

priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
               beta.normal = list(mean = 0, var = 2.72))

n.samples <- 5000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

out <- spPGOcc(occ.formula = ~ Community, 
             det.formula = ~ 1, 
             data = speciesList, 
             inits = inits, 
             priors = priors, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains,
             n.batch = 400, 
             batch.length = 25)
summary(out)

# Perform a posterior predictive check to assess model fit. 
ppc.out.sp <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out.sp)

# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
waicOcc(out)

# Concise summary of main parameter estimates
summary(out)
# Take a look at objects in resulting object
names(out)
str(out$beta.samples)
# Probability the effect of this variable on occupancy is positive
mean(out$beta.samples[, 2] > 0)
# Create simple plot summaries using MCMCvis package.

require(MCMCvis)
# Occupancy covariate effects ---------
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))


############## plotting predictions

river.pred <- seq(from = 0, 
                  to = max(siteCovariate$River)+2,
                  length.out = nrow(stations))

 community.pred <- seq(from = 0, 
                      to = max(siteCovariate$Community)+2,
                      length.out = nrow(stations))

# hunting.pred <- rep(1, nrow(stations))

trail.pred <- seq(from = min(siteCovariate$Trail.Distance), 
                  to = max(siteCovariate$Trail.Distance)+1,
                  length.out = nrow(stations))


####### PREDICTING ACROSS SPECTRUM OF DIST TO COMMUNITY VALUES ######
X.0 <- data.frame(
  intercept = 1,
  Community = community.pred)
X.0 <- as.matrix(X.0)

# out.pred <- predict(out, X.0, coords.0 = coords)
# 
# 
# plot.dat <- data.frame(x = stations$x, 
#                        y = stations$y, 
#                        mean.psi = apply(out.pred$psi.0.samples, 2, mean), 
#                        sd.psi = apply(out.pred$psi.0.samples, 2, sd), 
#                        stringsAsFactors = FALSE)

# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
library(stars)
library(ggplot2)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(title = 'Mean occurrence probability (community)') +
  theme_bw()


# begin interpolating values across landscape
coordsSP <- coords # copy
library(rgdal)
coordinates(coordsSP) <- c("x", "y")
proj4string(coordsSP) <- CRS("+init=epsg:24817")
coordsSP <- coordinates(coordsSP)

# make into a spdf
spdf <- SpatialPointsDataFrame(coordsSP, siteCovariate)
proj4string(spdf) <- CRS("+init=epsg:24817")
v <- terra::voronoi(spdf)

# nearest neighbor
spatVector <- vect(spdf, crs = "+init=epsg:24817")
ecuadorProj <- "+init=epsg:24817"
v <- voronoi(spatVector)
plot(v)
points(spatVector)

# crop it down
ext <- extent(c(xmin = min(coords$x)-0.05, 
                xmax = max(coords$x)+0.05, 
                ymin = max(coords$y)-0.15, 
                ymax = max(coords$y)+0.05))
croppedTiles <- crop(v, ext)
plot(croppedTiles, "River")
points(spatVector)
plot(croppedTiles, "Community")
points(spatVector)
# crs(croppedTiles)

# convert to raster
r <- rast(croppedTiles, res=10000)
print(r) # issue with these dimensions (should be more than 1 column and 1 row)
blankRaster <- terra::project(x = r, ecuadorProj)
rast <- terra::rasterize(croppedTiles, blankRaster)
plot(rast)
points(spatVector)

####### this is where I am stuck (found niche python article saying i need to increase significant figures, but unsure how to do that)
# basically my issue is trying to get this SpatVector into a raster with smaller tiles so I can use https://rspatial.org/analysis/4-interpolation.html to interpolate/krig




# tutorial for predicting unmarked objects: https://cran.rstudio.com/web/packages/unmarked/vignettes/spp-dist.html 


#################################### other attempts

# tutorial from https://rspatial.org/analysis/4-interpolation.html 

minimalCov <- data.frame(x = coords$x,
                         y = coords$y,
                         Community = siteCovariate$Community)
                         #River = siteCovariate$River)
occuPred <- predict(out,
                    X.0 = data.frame(Intercept = 1,
                                     Community = minimalCov$Community),
                    coords.0 = coords)
# occuPred <- predict(out,
#                     type = "state",
#                     newdata = minimalCov,
#                     na.rm = TRUE,
#                     inf.rm = TRUE) # https://doi90.github.io/lodestar/fitting-occupancy-models-with-unmarked.html#prediction 


levelplot(Predicted ~ minimalCov$x + minimalCov$y,
          data = occuPred,
          col.regions = rev(terrain.colors(100)),
          at = seq(0,1,length.out=101))








###################################################
##################### DEER ########################
###################################################

# load data
species <- read.csv("Deer.csv")
species <- data.frame(species[,-1], row.names=species[,1]) #rownames = stations

# occupancy 
library(unmarked)

# data as matrix
y = as.matrix(species)
y = y[ order(as.numeric(row.names(y))), ] #order matters

# occupancy covariates
siteCovariate = read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)

# unmarked df
ufo = unmarkedFrameOccu(y, 
                        siteCovs = siteCovaria te,
                        obsCovs = NULL)

plot(ufo)

# models (only Community significant)
Null = occu( ~1 ~1, ufo) # detection and occupancy are constant/unvarying

backTransform(Null, type = "det") # p AKA detection
backTransform(Null, type = "state") # psi AKA occupancy

Community = occu( ~1 ~CR, ufo) # detection is constant, but occupancy varies with the distance to the community
Community # only model that significantly fits the data
backTransform(Community, type = "det") # p AKA detection
# backTransform(Community, type = "state") # psi AKA occupancy 
# can't do because covariates are present -> would need to assign estimate for covariate(s)

River = occu( ~1 ~RR, ufo) # detection is constant, but occupancy varies with the distance to the river 
River

CommunityHabitat = occu( ~1 ~CR + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHabitat

CommunityHunted = occu( ~1 ~CR + Hunting, ufo) # occupancy varies with dist to community and hunting
CommunityHunted

# AIC values
BestModel = fitList(Null, Community, River, 
                    CommunityHabitat, CommunityHunted)
modSel(BestModel)

# predict
newdata = data.frame(0:10)
colnames(newdata)[1] = "RR"
predicted = predict(River, type="state", newdata=newdata, appendData=TRUE) # state = occupancy
par(pty="s")

# plot
library(ggplot2)
ggplot(NULL, aes(x=RR, y=Predicted)) + 
  geom_line(data=predicted, linetype="solid", size=1) +
  geom_ribbon(data=predicted, aes(ymin=lower, ymax=upper), fill="black", alpha=0.15) +
  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=1, expand = c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0:1), expand = c(0,0)) +
  theme_minimal() +
  coord_fixed(10) +
  labs(title = predicted$English, subtitle = predicted$Scientific, 
       x = "Distance from River (km)", 
       y = "Deer Occurrence Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


## Spatial occupancy model

stations <- read.csv("Stations2018.csv")
coords <- stations %>% 
  group_by(Station) %>%
  slice(1) %>%
  arrange(Station) %>%
  select(Station, y, x)

coords <- coords[,2:3]


library(spOccupancy)

speciesList <- list(y = species,
                    occ.covs = siteCovariate,
                    det.covs = siteCovariate,
                    coords = coords)

# Fit a spatial, single-species occupancy model
out.sp <- spPGOcc(occ.formula = ~RR, 
                  det.formula = ~ 1, 	          
                  data = speciesList, 
                  n.batch = 400, 
                  batch.length = 25,
                  n.thin = 5, 
                  n.burn = 5000, 
                  n.chains = 3,
                  n.report = 100)
summary(out.sp)




# Format with abbreviated specification of inits for alpha and beta.
inits <- list(alpha = 0, 
              beta = 0, 
              z = apply(speciesList$y, 1, max, na.rm = TRUE))

priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
               beta.normal = list(mean = 0, var = 2.72))

n.samples <- 5000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

out <- PGOcc(occ.formula = ~ RR, 
             det.formula = ~ 1, 
             data = speciesList, 
             inits = inits, 
             n.samples = n.samples, 
             priors = priors, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)


# Perform a posterior predictive check to assess model fit. 
ppc.out.sp <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out.sp)

# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
waicOcc(out)

# Concise summary of main parameter estimates
summary(out)
# Take a look at objects in resulting object
names(out)
str(out$beta.samples)
# Probability the effect of this variable on occupancy is positive
mean(out$beta.samples[, 2] > 0)
# Create simple plot summaries using MCMCvis package.

require(MCMCvis)
# Occupancy covariate effects ---------
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))


############## plotting predictions

river.pred <- seq(from = min(siteCovariate$River)-100, 
                  to = max(siteCovariate$River)+100,
                  length.out = nrow(stations))

community.pred <- seq(from = min(siteCovariate$Community)-100, 
                      to = max(siteCovariate$Community)+100,
                      length.out = nrow(stations))

# hunting.pred <- rep(1, nrow(stations))

trail.pred <- seq(from = min(siteCovariate$Trail.Distance), 
                  to = max(siteCovariate$Trail.Distance)+1,
                  length.out = nrow(stations))


####### PREDICTING ACROSS SPECTRUM OF DIST TO RIVER VALUES ######
X.0 <- data.frame(
  intercept = 1,
  River = river.pred)
X.0 <- as.matrix(X.0)

out.pred <- predict(out, X.0)


plot.dat <- data.frame(x = stations$x, 
                       y = stations$y, 
                       mean.psi = apply(out.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.pred$psi.0.samples, 2, sd), 
                       stringsAsFactors = FALSE)

# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
library(stars)
library(ggplot2)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(title = 'Mean occurrence probability (rivers)') +
  theme_bw()














###################################################
##################### PACA ########################
###################################################

# load data
species <- read.csv("Paca.csv")
species <- data.frame(species[,-1], row.names=species[,1]) #rownames = stations

# occupancy 
library(unmarked)

# data as matrix
y = as.matrix(species)
y = y[ order(as.numeric(row.names(y))), ] #order matters

# occupancy covariates
siteCovariate = read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)


# unmarked df
ufo = unmarkedFrameOccu(y, 
                        siteCovs = siteCovariate,
                        obsCovs = NULL)

plot(ufo)

# models (non significant)
Null = occu( ~1 ~1, ufo) # detection and occupancy are constant/unvarying

Community = occu( ~1 ~CR, ufo) # detection is constant, but occupancy varies with the distance to the community
Community 
backTransform(Community, type = "det") # p AKA detection
# backTransform(Community, type = "state") # psi AKA occupancy 
# can't do because covariates are present -> would need to assign estimate for covariate(s)

River = occu( ~1 ~RR, ufo) # detection is constant, but occupancy varies with the distance to the river 
River

CommunityHabitat = occu( ~1 ~CR + Habitat, ufo) # occupancy varies with dist to community and habitat
CommunityHabitat

CommunityHunted = occu( ~1 ~CR + Hunting, ufo) # occupancy varies with dist to community and hunting
CommunityHunted

# AIC values
BestModel = fitList(Null, Community, River, 
                    CommunityHabitat, CommunityHunted)
modSel(BestModel)

# predict
newdata = data.frame(0:10)
colnames(newdata)[1] = "RR"
predicted <- predict(River, type="state", newdata=newdata, appendData=TRUE) # state = occupancy
par(pty="s")

# plot
library(ggplot2)
ggplot(NULL, aes(x=RR, y=Predicted)) + 
  geom_line(data=predicted, linetype="solid", size=1) +
  geom_ribbon(data=predicted, aes(ymin=lower, ymax=upper), fill="black", alpha=0.15) +
  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=1, expand = c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0:1), expand = c(0,0)) +
  theme_minimal() +
  coord_fixed(10) +
  labs(title = predicted$English, subtitle = predicted$Scientific, 
       x = "Distance from River (km)", 
       y = "Paca Occurrence Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


## Spatial occupancy model

stations <- read.csv("Stations2018.csv")
coords <- stations %>% 
  group_by(Station) %>%
  slice(1) %>%
  arrange(Station) %>%
  select(Station, y, x)

coords <- coords[,2:3]


library(spOccupancy)

speciesList <- list(y = species,
                       occ.covs = siteCovariate,
                       det.covs = siteCovariate,
                       coords = coords)

# Fit a spatial, single-species occupancy model
out.sp <- spPGOcc(occ.formula = ~CR + Habitat, 
                  det.formula = ~ 1, 	          
                  data = speciesList, 
                  n.batch = 400, 
                  batch.length = 25,
                  n.thin = 5, 
                  n.burn = 5000, 
                  n.chains = 3,
                  n.report = 100)
summary(out.sp)




# Format with abbreviated specification of inits for alpha and beta.
inits <- list(alpha = 0, 
                   beta = 0, 
                   z = apply(speciesList$y, 1, max, na.rm = TRUE))

priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal = list(mean = 0, var = 2.72))

n.samples <- 5000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

out <- PGOcc(occ.formula = ~ CR + Habitat, 
             det.formula = ~ 1, 
             data = speciesList, 
             inits = inits, 
             n.samples = n.samples, 
             priors = priors, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)


# Perform a posterior predictive check to assess model fit. 
ppc.out.sp <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out.sp)

# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
waicOcc(out)

# Concise summary of main parameter estimates
summary(out)
# Take a look at objects in resulting object
names(out)
str(out$beta.samples)
# Probability the effect of this variable on occupancy is positive
mean(out$beta.samples[, 2] > 0)
# Create simple plot summaries using MCMCvis package.

require(MCMCvis)
# Occupancy covariate effects ---------
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))


############## plotting predictions
str(hbefElev)
str(siteCovariate)
range(siteCovariate$River)

river.pred <- seq(from = min(siteCovariate$River)-100, 
                  to = max(siteCovariate$River)+100,
                  length.out = nrow(stations))

community.pred <- seq(from = min(siteCovariate$Community)-100, 
                         to = max(siteCovariate$Community)+100,
                         length.out = nrow(stations))

# hunting.pred <- rep(1, nrow(stations))

trail.pred <- seq(from = min(siteCovariate$Trail.Distance), 
                      to = max(siteCovariate$Trail.Distance)+1,
                      length.out = nrow(stations))


####### PREDICTING ACROSS SPECTRUM OF DIST TO COMMUNITY VALUES ######
X.0 <- data.frame(
  intercept = 1,
  Community = community.pred,
  HabitatSwamp = 0,
  HabitatUpland = 0
)
X.0 <- as.matrix(X.0)
  
out.pred <- predict(out, X.0)


plot.dat <- data.frame(x = stations$x, 
                       y = stations$y, 
                       mean.psi = apply(out.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.pred$psi.0.samples, 2, sd), 
                       stringsAsFactors = FALSE)

# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
library(stars)
library(ggplot2)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(title = 'Mean occurrence probability (communities)') +
  theme_bw()



####### PREDICTING ACROSS SPECTRUM OF DIST TO COMMUNITY VALUES ######
X.0 <- data.frame(
  intercept = 1,
  River = mean(siteCovariate$River),
  Community = community.pred,
  Hunting = 1,
  HabitatSwamp = 0,
  HabitatUpland = 0,
  Trail.Distance = mean(siteCovariate$Trail.Distance)
)
X.0 <- as.matrix(X.0)

out.pred <- predict(out, X.0)


plot.dat <- data.frame(x = stations$x, 
                       y = stations$y, 
                       mean.psi = apply(out.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.pred$psi.0.samples, 2, sd), 
                       stringsAsFactors = FALSE)

# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
library(stars)
library(ggplot2)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(title = 'Mean occurrence probability (rivers)') +
  theme_bw()












