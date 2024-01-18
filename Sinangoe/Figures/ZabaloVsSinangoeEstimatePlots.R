# Zabalo vs Sinangoe occupancy estimates plot

# load packages
require(ggplot2)

# load data
load('occupancyEstimates.RData')

# work around to make x-axis labels with italics and divided into two lines
peccary <- ~ atop(paste("Collared peccary"), paste("(", italic("Dicotyles tajacu"), ")"))
paca <- ~ atop(paste("Lowland paca"), paste("(", italic("Cuniculus paca"), ")"))
trumpeter <- ~ atop(paste("Grey-winged trumpeter"), paste("(", italic("Psophia crepitans"), ")"))

# plot
ggplot(occupancyEstimates, aes(x = factor(Species, levels = unique(Species)),
                               y = Predicted,
                               fill = factor(Community, levels = unique(Community)))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("darkorange", "deepskyblue")) +
  labs(x = "Species", y = "Average occupancy probability") +
  theme_classic() +
  theme(text = element_text(family = "Times", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c(peccary, paca, trumpeter))

# save it
ggsave(filename = "ZabaloSinangoeOccupancyEstimates.png", 
       width = 8, height = 4)

