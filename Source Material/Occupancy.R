# basica occupancy model

# load data
species = read.csv("CollaredPeccary.csv")
species = data.frame(species[,-1], row.names=species[,1]) #rownames = stations

# cccupancy 
library(unmarked)

# data as matrix
y = as.matrix(species)
y = y[ order(as.numeric(row.names(y))), ] #order matters

# occupancy covariates
siteCovariate = read.csv("siteCovs2018.csv")

# unmarked df
ufo = unmarkedFrameOccu(y, 
                        siteCovs = siteCovariate,
                        obsCovs = NULL)

plot(ufo)

# models
Null = occu( ~1 ~1, ufo)
Community = occu( ~1 ~CR, ufo)
River = occu( ~1 ~RR, ufo)

CommunityHabitat = occu( ~1 ~CR + Habitat, ufo)
CommunityHunted = occu( ~1 ~CR + Hunting, ufo)

# AIC values
BestModel = fitList(Null, Community, River, 
                    CommunityHabitat, CommunityHunted)
modSel(BestModel)

# predict
newdata = data.frame(0:10)
colnames(newdata)[1] = "CR"
predicted = predict(Community, type="state", newdata=newdata, appendData=TRUE) # state = occupancy
par(pty="s")

# plot
library(ggplot2)
ggplot(NULL, aes(x=CR, y=Predicted)) + 
  geom_line(data=predicted, linetype="solid", size=1) +
  geom_ribbon(data=predicted, aes(ymin=lower, ymax=upper), fill="black", alpha=0.15) +
  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=1, expand = c(0,0)) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0:1), expand = c(0,0)) +
  theme_minimal() +
  coord_fixed(10) +
  labs(title = predicted$English, subtitle = predicted$Scientific, 
       x = "Distance from Community (km)", 
       y = "Occurrence Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
