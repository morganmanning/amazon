# conducting a power analysis based on this Ken Kellner's tutorial: 
# https://cran.r-project.org/web/packages/unmarked/vignettes/powerAnalysis.html 

# setwd
setwd("~/Documents/amazon/Global/Data")

# packages
require(lubridate)
require(unmarked)
require(camtrapR)
require(kableExtra)

# load in all necessary data
Data <- read.csv("AllIndependentRecordsFormatted.csv")
Data$DateTimeOriginal <- parse_date_time(Data$DateTimeOriginal, c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S"))
Traps <- read.csv("AllStationsFormatted.csv")
siteCovariate <- read.csv("AllCommunityCovariates.csv")
species <- c("Pecari tajacu", "Mazama americana", "Cuniculus paca", "Psophia crepitans")

# camera operability matrix
Operation <- cameraOperation(CTtable = Traps,
                             stationCol = "Station",
                             cameraCol = "Camera",
                             setupCol = "Setup_date",
                             retrievalCol = "Retrieval_date",
                             hasProblems = TRUE,
                             byCamera = FALSE,
                             allCamsOn = FALSE,
                             camerasIndependent = FALSE,
                             dateFormat = "%Y-%m-%d",
                             writecsv = FALSE)
occasion <- 10

for (i in 1:length(species)) {
DetHis = detectionHistory(recordTable = Data,
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
                          species = species[i]) #change species here
justDetHis <- DetHis[["detection_history"]]

# into correct format
ufo <- unmarkedFrameOccu(justDetHis,
                         siteCovs = siteCovariate,
                         obsCovs = NULL)

# global model
output <- occu(~Community ~Community + Temperature + percentNatural, 
               ufo, 
               control = 10000,
               starts = c(rep(0, times = 10)))
powerAnalysis(output)

}





output <- occu(~1 ~1, ufo, 
               control = 10000,
               starts = c(0, 0)) # starting values for parameters


powerAnalysis(output)
effect_sizes <- list(state = c(intercept = 0, CommunitySinangoe = , CommunitySiona =, CommunityZabalo = , Temperature = , percentNatural = ), 
                     det=c(intercept=0, CommunitySinangoe = , CommunitySiona = , CommunityZabalo = ))
(pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05, nsim=20))


# 50 sites and 3 obs per site
(pa2 <- powerAnalysis(template_model, effect_sizes, design=list(M=50, J=3), nsim=20))