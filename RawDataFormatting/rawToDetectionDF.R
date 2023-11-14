rm(list = ls())

############################################################################
###################### FORMATTING DETECTION MATRICES #######################
############################################################################
# requirements: 
# notes:


setwd("~/Documents/amazon/RawDataFormatting")

library(openxlsx)
library(camtrapR)
library(dplyr)
require(unmarked)



############################################################################
############################# LOAD DATA ####################################
############################################################################

# Sinangoe
Data <- read.csv("Sinangoe/RecordTable.csv") # this for camtrapR (Steps 1-5)
Traps <- read.csv("Sinangoe/Stations.csv")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---------------------- CAMERA INDEPENDENCE ----------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# create DateTimeOrginal column in proper format
Data$DateTimeOriginal = strptime(paste(as.Date(Data$Date, format = "%d/%m/%Y"), # changed this from OG code
                                       Data$Time),
                                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# temporal independence
source("assessTemporalIndependence.R")
Time30 = assessTemporalIndependence(intable = Data,
                                    deltaTimeComparedTo = "lastIndependentRecord",
                                    columnOfInterest = "Species",
                                    stationCol = "Station",
                                    cameraCol = "CameraName",
                                    camerasIndependent = FALSE,
                                    minDeltaTime = 30)

# check that the previous function worked (e.g. no duplicates)
Time30 = Time30[order(Time30$DateTimeOriginal),]
Time30 <- Time30 %>% 
  distinct(.keep_all = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ----------------------------- DETECTIONS ------------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# explore total detections
total = data.frame(table(Time30$Species)) # number of detections / species
colnames(total) = c("Species", "Total")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# -------------------------- OCCUPANCY SET-UP ---------------------------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# camera operability matrix
Operation = cameraOperation(CTtable = Traps,
                            stationCol = "Station",
                            cameraCol = "Camera",
                            setupCol = "Setup_date",
                            retrievalCol = "Retrieval_date",
                            hasProblems = TRUE,
                            byCamera = FALSE,
                            allCamsOn = FALSE,
                            camerasIndependent = FALSE,
                            dateFormat = "%d/%m/%Y",
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
