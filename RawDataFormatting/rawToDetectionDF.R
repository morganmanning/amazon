rm(list = ls())

setwd("/Users/mesbach/Dropbox (Princeton)/2. Publications/Camera Traps")

year=2018

# load data
library(openxlsx)
Data = read.csv(paste0("Data/RecordTable",year,".csv")) # this for camtrapR (Steps 1-5)
Traps = read.csv(paste0("Data/Stations",year,".csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---------------------- CAMERA INDEPENDENCE ----------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(camtrapR)

# create DateTimeOrginal column in proper format
Data$DateTimeOriginal = strptime(paste(as.Date(Data$Date, format = "%m/%d/%y"),
                                       Data$Time),
                                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# temporal independence
source("Code/Functions/assessTemporalIndependence.R")
Time30 = assessTemporalIndependence(intable = Data,
                                    deltaTimeComparedTo = "lastIndependentRecord",
                                    columnOfInterest = "Species",
                                    stationCol = "Station",
                                    cameraCol = "CameraName",
                                    camerasIndependent = FALSE,
                                    minDeltaTime = 30)

# check that the previous function worked (e.g. no duplicates)
Time30 = Time30[order(Time30$DateTimeOriginal),]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ----------------------------- DETECTIONS ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# explore total detections
total = data.frame(table(Time30$Species)) # number of detections / species
colnames(total) = c("Species", "Total")

Hunting = subset(Time30, Station %in% c(1, 2, 5, 10, 11, 14, 15, 16, 19, 20, 23, 24, 25, 27, 28))
NonHunting = subset(Time30, Station %in% c(3, 4, 6, 7, 8, 9, 12, 13, 17, 18, 21, 22, 26, 29, 30))

# Hunting detections
hunted = data.frame(table(Hunting$Species)) # number of detections / species
colnames(hunted) = c("Species", "Hunted")

# NonHunting detections
sepicho = data.frame(table(NonHunting$Species)) # number of detections / species
colnames(sepicho) = c("Species", "Sepicho")

detections = merge(total, hunted, by="Species", all=TRUE)
detections = merge(detections, sepicho, by="Species", add=TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ------------------------- ACTIVITY PATTERNS ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# put time back into character format
Time30$DateTimeOriginal = as.character(Time30$DateTimeOriginal)

# single-species diel activity kernel density estimation plots
activityDensity(recordTable = Time30,
                species = "Pecari-tajacu",
                allSpecies = FALSE,
                speciesCol = "Species",
                recordDateTimeCol = "DateTimeOriginal",
                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                plotR = TRUE,
                writePNG = FALSE,
                plotDirectory,
                createDir = FALSE,
                pngMaxPix = 1000,
                add.rug = TRUE) 

# two-species diel activity overlap plots and estimates
activityOverlap(recordTable = Time30,
                speciesA = "Pecari-tajacu",
                speciesB = "Tayassu-pecari",
                speciesCol = "Species",
                recordDateTimeCol = "DateTimeOriginal",
                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                plotR = TRUE,
                writePNG = FALSE,
                addLegend = FALSE,
                legendPosition = "topleft",
                plotDirectory,
                createDir = FALSE,
                pngMaxPix = 1000,
                add.rug = TRUE,
                overlapEstimator = c("Dhat1", "Dhat4", "Dhat5"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ----------------------------- MAPPING ---------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(rgdal)

# backgroundPoylogon (for function below)
StudyArea = readOGR("/Users/mesbach/Dropbox (Princeton)/6. PhD/1. Dissertation/Part 3/4. Camera Traps/yr2018/1. GIS/StudyArea/Background.shp")

# species richness
detectionMaps(CTtable = Traps,
              recordTable = Time30,
              Xcol = "x", 
              Ycol = "y",
              backgroundPolygon = StudyArea,
              stationCol = "Station",
              speciesCol = "Species",
              speciesPlots = FALSE, #invividual species
              richnessPlot = TRUE, #richness across entire project
              printLabels = FALSE, #station numbers in red
              addLegend = FALSE) 

# species
unique(Time30$Species)
species = "Mazama-gouazoubira"

# species presence by station
detectionMaps(CTtable = Traps,
              recordTable = Time30,
              Xcol = "x", 
              Ycol = "y",
              backgroundPolygon = StudyArea,
              stationCol = "Station",
              speciesCol = "Species",
              speciesToShow = species,
              speciesPlots = TRUE, #invividual species
              richnessPlot = FALSE, #richness across entire project
              printLabels = FALSE,
              addLegend = FALSE) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- OCCUPANCY SET-UP ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# camera operability matrix
Operation = cameraOperation(CTtable = Traps,
                            stationCol = "Station",
                            cameraCol = "CameraName",
                            setupCol = "Setup_date",
                            retrievalCol = "Retrieval_date",
                            hasProblems = TRUE,
                            byCamera = FALSE,
                            allCamsOn = FALSE,
                            camerasIndependent = FALSE,
                            dateFormat = "%m/%d/%Y",
                            writecsv = FALSE)

# occasion length
occasion = 2

# species detection histories for occupancy analyses
DetHis = detectionHistory(recordTable = Time30,
                          camOp = Operation,
                          output = "binary", # binary or count
                          stationCol = "Station",
                          speciesCol = "Species",
                          recordDateTimeCol = "DateTimeOriginal",
                          recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                          day1 = "Station",
                          occasionLength = occasion,
                          datesAsOccasionNames = FALSE,
                          timeZone = "America/Guayaquil",
                          includeEffort = TRUE,
                          scaleEffort = FALSE,
                          #maxNumberDays = 90, #need to think about this
                          species = species) #change species here

# Deer = DetHis[["detection_history"]]
# write.csv(Deer, "~/Desktop/Occupancy/Deer.csv", row.names=T)


















#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- OCCUPANCY ANALYSIS -------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(unmarked)

y = as.matrix(DetHis$detection_history)
y = y[ order(as.numeric(row.names(y))), ]

# occupancy covariates
siteCovariate = read.csv(paste0("Data/siteCovs",year,".csv"))

# unmarked df
ufo = unmarkedFrameOccu(y, 
                        siteCovs = siteCovariate,
                        obsCovs = NULL)

plot(ufo)
summary(ufo)

# ----------------------- MODELS ----------------------------------------------#

# null
Null = occu( ~1 ~1, ufo)

# occupancy
Community = occu( ~1 ~CR, ufo)
River = occu( ~1 ~RR, ufo)

# multiple
CommunityHabitat = occu( ~1 ~CR + Habitat, ufo)
CommunityHunted = occu( ~1 ~CR + Hunting, ufo)

# detection
CommunityEffort = occu( ~ScaleEffort ~CR, ufo)

# AIC values
BestModel = fitList(Null, Community, River, 
                    CommunityHabitat, CommunityHunted, CommunityEffort)
modSel(BestModel)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- PLOT PREDICITIONS --------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# ONE VARIABLE
newdata = data.frame(0:10)
colnames(newdata)[1] = "CR"
predicted = predict(Community, type="state", newdata=newdata, appendData=TRUE) # state = occupancy
par(pty="s")

plot(Predicted~CR, predicted, type="l", ylim=0:1, lwd=2, 
     xlab = "Distance from Community (km)",
     ylab = "Occurance Probability", 
     main = paste0(species))
lines(upper~CR, predicted, lty=2, col="red")
# lines(lower~CR, predicted, lty=2, col="red")
# lines(newdata$CR, predicted[,"Predicted"]+1.96*predicted[,"SE"], lty=2, col="red")
lines(newdata$CR, predicted[,"Predicted"]-1.96*predicted[,"SE"], lty=2, col="red")

# ------------------------------ GGPLOT ---------------------------------------#

library(ggplot2)

# add restriction
predicted$Restriction = "Hunted"
predicted[6:11,6] = "Not Hunted"

# ensure proper order
predicted$Restriction = factor(predicted$Restriction, 
                               levels = c("Not Hunted","Hunted"))

# gap
predicted = rbind(predicted, predicted[rep(6, 1), ])
predicted[12, "Restriction"] = "Hunted"

# LCL
predicted$lower = predicted$Predicted - predicted$SE*1.96
predicted$lower[predicted$lower<0] = 0

# names
English = read.xlsx("/Users/mesbach/Dropbox (Princeton)/2. Publications/Camera Traps/Code/Names.xlsx")
predicted$English = English$English[match(species, English$Scientific)]
predicted$Scientific = species
predicted$Scientific = gsub("-"," ",as.character(predicted$Scientific))

# plot
ggplot(NULL, aes(x=CR, y=Predicted)) + 
  geom_line(data=predicted, linetype="solid", size=1, aes(color=Restriction)) +
  geom_ribbon(data=predicted, aes(ymin=lower, ymax=upper), fill="black", alpha=0.15) +
  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=1, expand = c(0,0)) + 
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0:1), expand = c(0,0)) + 
  scale_color_manual(name="Restriction", values=c("#009E73","#D55E00")) + 
  theme_minimal() +
  coord_fixed(10) +
  labs(title = predicted$English, subtitle = predicted$Scientific, 
       x = "Distance from Community (km)", 
       y = "Occurrence Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_vline(aes(xintercept = 5, linetype = "~5 km")) +
  scale_linetype_manual(name = "Hunting Limit", values = 5) + 
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))
