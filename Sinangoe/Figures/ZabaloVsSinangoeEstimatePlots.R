# Zabalo vs Sinangoe occupancy estimates plot

# load packages
require(ggplot2)
require(dplyr)
library(magick)
library(grid)

# load data
load('occupancyEstimates.RData')

# prep for dot and whisker plot
meanSE <- occupancyEstimates |>
  group_by(Community, Species) |>
  summarize(mean = mean(Predicted), se = sd(Predicted)/sqrt(n()), sd = sd(Predicted))
meanSE$Species <- factor(meanSE$Species, 
                         levels=c("Collared peccary", "Lowland paca", "Grey-winged trumpeter"))

# work around to make x-axis labels with italics and divided into two lines
peccary <- ~ atop(paste("Collared peccary"), paste("(", italic("Dicotyles tajacu"), ")"))
paca <- ~ atop(paste("Lowland paca"), paste("(", italic("Cuniculus paca"), ")"))
trumpeter <- ~ atop(paste("Grey-winged trumpeter"), paste("(", italic("Psophia crepitans"), ")"))

# plot it
dodge <- position_dodge(width = 0.25)
plot <- ggplot(meanSE, aes(x = Species,
                   y = mean,
                   color = Community)) +
  geom_point(aes(color = Community), position = dodge) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = Community), 
                position = dodge, width = 0.1) +
  scale_color_manual(values = c("darkorange", "royalblue")) +
  scale_x_discrete(labels = c(peccary, paca, trumpeter)) +
  labs(x = "Species", y = "Probability of Occupancy") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position="top")


f2 <- image_read("peccary.png")
f2 <- image_trim(f2)
f3 <- image_read("paca.png")
f3 <- image_trim(f3)
f4 <- image_read("trumpeter.png")
f4 <- image_trim(f4)
plot
grid.raster(f2,x=0.22,y=0.14, width=0.1)
grid.raster(f3,x=0.52,y=0.14, width=0.1)
grid.raster(f4,x=0.81,y=0.20, width=0.1)

# save it
ggsave(filename = "ZabaloSinangoeOccupancyEstimates.tiff", 
       width = 8, height = 4)

