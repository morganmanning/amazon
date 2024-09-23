# Zabalo vs Sinangoe occupancy estimates plot

# load packages
require(ggplot2)
require(dplyr)
require(rphylopic)

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


# animal silhouettes 
peccPic <- get_phylopic(uuid = get_uuid(name = "Dicotyles tajacu", n = 1))
pacaPic <- get_phylopic(uuid = get_uuid(name = "Cuniculus paca", n = 1))
trumpPic <- get_phylopic(uuid = get_uuid(name = "Psophia crepitans", n = 1))

# plot it
dodge <- position_dodge(width = 0.25)
transparency <- 0.2
ggplot(meanSE, aes(x = Species,
                   y = mean,
                   color = Community)) +
  geom_point(aes(color = Community, shape = Community), position = dodge, size = 1.75) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = Community, linetype = Community),
                position = dodge, width = 0.1) +
  #scale_color_manual(values = c("darkorange", "royalblue")) +
  scale_color_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(16, 15)) +
  scale_x_discrete(labels = c(peccary, paca, trumpeter)) +
  labs(x = "Species", y = "Probability of occupancy") +
  ylim(c(0,1)) +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position="top") +
  add_phylopic(peccPic, alpha = transparency, x = 1, y = 0.07, ysize = 0.2) +
  add_phylopic(pacaPic, alpha = transparency, x = 2, y = 0.07, ysize = 0.175) +
  add_phylopic(trumpPic, alpha = transparency, x = 3, y = 0.1, ysize = 0.275)


  

# save it
ggsave(filename = "ZabaloSinangoeOccupancyEstimates.png", 
       width = 8, height = 4)

