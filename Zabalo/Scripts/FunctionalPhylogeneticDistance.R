
##### CALCULATE FUNCTIONAL DIVERSITY #####
setwd("/Users/morganmanning/Documents/amazon/Zabalo/Data")
# load data
peccary <- read.csv("CollaredPeccary.csv")
peccary <- data.frame(peccary[,-1], row.names=peccary[,1]) #rownames = stations

deer <- read.csv("Deer.csv")
deer <- data.frame(deer[,-1], row.names=deer[,1]) #rownames = stations

paca <- read.csv("Paca.csv")
paca <- data.frame(paca[,-1], row.names=paca[,1]) #rownames = stations

# data as matrix
peccaryOccupancy <- as.matrix(peccary)
peccaryOccupancy <- peccaryOccupancy[ order(as.numeric(row.names(peccaryOccupancy))), ] #order matters

deerOccupancy <- as.matrix(deer)
deerOccupancy <- deerOccupancy[ order(as.numeric(row.names(deerOccupancy))), ] #order matters

pacaOccupancy <- as.matrix(paca)
pacaOccupancy <- pacaOccupancy[ order(as.numeric(row.names(pacaOccupancy))), ] #order matters


# occupancy covariates
siteCovariate = read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$X)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$Community <- siteCovariate$Community/1000
siteCovariate$River <- siteCovariate$River/1000

# Read in the animal-specific data
masterlist<-read.csv("./FPDistData/master_species_list_updated_7April2014.csv",h=T) #master list from Beaudrot
pantheria <- read.delim("./FPDistData/PanTHERIA_1-0_WR05_Aug2008.txt", header = TRUE, sep = "\t") # https://esapubs.org/archive/ecol/E090/184/metadata.htm

for ( col in 1:ncol(pantheria)){
  colnames(pantheria)[col] <-  gsub("^.*?\\_", "", colnames(pantheria)[col])
} # removing wacky column beginnings

# replace all -999 with NA
pantheria[pantheria == -999] <- NA 


pantheria$Binomial <- gsub(" ", "_", pantheria$Binomial)
rownames(pantheria) <- pantheria$Binomial


# pull the desired traits
require(FD)
require(dplyr)
# names(pantheria)
traits <- pantheria %>%
  select(AdultHeadBodyLen_mm,
         AdultBodyMass_g, 
         LitterSize, 
         GR_Area_km2, # geographic range area
         ActivityCycle, # (1) nocturnal only, (2) nocturnal/crepuscular, cathemeral, crepuscular or diurnal/crepuscular and (3) diurnal only
         TrophicLevel, # (1) herbivore (not vertebrate and/or invertebrate), (2) omnivore (vertebrate and/or invertebrate plus any of the other categories) and (3) carnivore (vertebrate and/or invertebrate only)
         HabitatBreadth, # number of habitats used: above ground dwelling, aquatic, fossorial and ground dwelling
         DietBreadth, # number of diets used: vertebrate, invertebrate, fruit, flowers/nectar/pollen, leaves/branches/bark, seeds, grass and roots/tubers
         Order,
         Family,
         Genus,
         Species)

# reclassify factors
traits$TrophicLevel <- ifelse(traits$TrophicLevel == 1, 'herbivore', 
                              ifelse(traits$TrophicLevel == 2, 'omnivore', 'carnivore'))
traits$ActivityCycle <- ifelse(traits$ActivityCycle == 1, 'nocturnal', 
                              ifelse(traits$ActivityCycle == 2, 'other', 'diurnal'))
traits$HabitatBreadth <- as.factor(traits$HabitatBreadth)
traits$DietBreadth <- as.factor(traits$DietBreadth)
traits$ActivityCycle <- as.factor(traits$ActivityCycle)
traits$TrophicLevel <- as.factor(traits$TrophicLevel)
str(traits)

# check out the animals of interest
(check <- traits[rownames(traits) == 'Pecari_tajacu' | 
                   rownames(traits) == 'Mazama_gouazoubira' | 
                   rownames(traits) == 'Cuniculus_paca' ,])

# calculate functional distance
FDist <- as.matrix(gowdis(traits))


FDist[,c('Pecari_tajacu','Mazama_gouazoubira', 'Cuniculus_paca')]

# export large matrix of similarity
save(FDist, file = 'FuntionalDistance.RData')




# look at the nearest five animals
nSimilar <- 5

# pull the most similar 5 animals for each animal and put it in a df
mostSimilar <- data.frame()
for (i in 1:ncol(FDist)){
  speciesOfInterest <- data.frame(species = colnames(FDist)[i],
                                  nearest = rownames(FDist),
                                  functionalDist = FDist[,i])
  speciesOfInterest <- speciesOfInterest[order(speciesOfInterest$functionalDist, 
                                               decreasing = FALSE),] 
  speciesOfInterestTop <- speciesOfInterest[2:(nSimilar+1),] # don't pull comparison of same spp.
  mostSimilar <- rbind(mostSimilar, speciesOfInterestTop)
}
rownames(mostSimilar) <- 1:nrow(mostSimilar)


























# # Only keep the animals of interest
# # Collared peccary: Pecari tajacu
# # Brown brocket: Mazama gouazoubira
# # Paca: Cuniculus
# studyAnimals <- pantheria[(pantheria$Genus == 'Pecari' & pantheria$Species == 'tajacu') |
#                             (pantheria$Genus == 'Cervus' & pantheria$Species == 'elaphus') |
#                             (pantheria$Genus == 'Cuniculus' & pantheria$Species == 'paca'),]
# 
# studyAnimals$Guild <- c("Omnivore", "Herbivore", "Herbivore")
# studyAnimals$Class <- "Mammalia"
# # head(studyAnimals)
# 
# 
# # using only traits that Lydia used
# spTraits <- data.frame(# Binomial = studyAnimals$Binomial, 
#                   # Class = studyAnimals$Class,
#                   # Family = studyAnimals$Family,
#                   BodyLength = studyAnimals$AdultHeadBodyLen_mm, 
#                   LitterSize = studyAnimals$LitterSize, 
#                   GR_Area = studyAnimals$GR_Area_km2, 
#                   ActivityCycle = as.factor(studyAnimals$ActivityCycle), 
#                   HabitatBreadth = as.factor(studyAnimals$HabitatBreadth), 
#                   DietBreadth = as.factor(studyAnimals$DietBreadth), 
#                   Guild = as.factor(studyAnimals$Guild))
# rownames(spTraits) <- c("Deer", "Paca", "Peccary")
# 
# # add category for each variable
# str(spTraits)
# traitCategory <- data.frame(trait_name = colnames(spTraits),
#                             trait_type = c("Q", "Q", "Q", "N", "N", "N", "N"))
# 
# 
# # species presence at each site
# # rows as sites, species as columns with names matching spTraits row names
# sitePresence <- data.frame(Site = 1:30,
#                            Deer = NA,
#                            Peccary = NA,
#                            Paca = NA)
# for (i in 1:nrow(peccaryOccupancy)) {
#   sitePresence$Peccary[i] <- ifelse(sum(peccaryOccupancy[i,], na.rm = TRUE) == 0, 0, 1)
# }
# 
# for (i in 1:nrow(pacaOccupancy)) {
#   sitePresence$Paca[i] <- ifelse(sum(pacaOccupancy[i,], na.rm = TRUE) == 0, 0, 1)
# }
# 
# for (i in 1:nrow(deerOccupancy)) {
#   sitePresence$Deer[i] <- ifelse(sum(deerOccupancy[i,], na.rm = TRUE) == 0, 0, 1)
# }
# sitePresence <- sitePresence[,-1]
# sitePresence <- sitePresence[rowSums(sitePresence[])>0,] # removing all sites where no species were found
# 
# spDist <- mFD::funct.dist(sp_tr = spTraits,
#                           tr_cat = traitCategory,
#                           metric = "gower") # since I have categorical variables
# spDist
                          

################ GENETIC DISTANCE ###############
# cophenetic.phylo() computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
# phylogenetic tree site: http://vertlife.org/phylosubsets/ 
  # tree-pruner-1b0220f3-d821-4782-ae78-d1403c4070dd

# mammal tree: https://github.com/lbeaudrot/Elevational-Shifts/blob/master/mammalST_MSW05_all.tre 
# https://raw.githubusercontent.com/lbeaudrot/Elevational-Shifts/master/FritzTree.rs200k.100trees.tre 

# Create phylogenetic and trail matrices to use as input for FPDist

library(picante)
library(ecodist)
library(vegan)
library(FD)

#input data
# mat <- read.delim("./FPDistData/mammal_com.csv", sep=",", row.names=1) # doesn't include red deer
com <- read.delim("./FPDistData/species_list.csv", sep="", row.names=1) # Latin names with _ as space
com <- data.frame(sp = c('Pecari_tajacu', 'Mazama_gouazoubira', 'Cuniculus_paca'))
# traits <- read.delim("mammal_traits.csv", sep=",", row.names=1) # doesn't include red deer

#load unresolved phylo
phylo <- read.nexus("./FPDistData/mammalST_MSW05_all.tre")
phylo <- phylo[[1]]

PDist <- cophenetic.phylo(phylo)


# 
# 
# #load resolved phylo
phylo_res <- read.nexus("./FPDistData/FritzTree.rs200k.100trees.tre")
phylo_names <- read.csv("./FPDistData/global_tip_labels_mod.csv")
phylo_names1 <- as.character(phylo_names$sp)
# 
# ### match phylo data with incidence
# phylomat<-list()
# for (i in 1:nrow(phylo_res)) {
#   my_phylo_res<-phylo_res[[i]]
#   #modify species names in phylo to match TEAM data
#   #write.table(phylo_res$tip.label,file="globa_tip_labels.txt")
#   my_phylo_res$tip.label<-phylo_names1
#   
#   ### prune phylo to match camera trap data
#   myphylosp <- match.phylo.data(my_phylo_res, com)
#   phy <- myphylosp$phy
#   phydist <- as.data.frame(cophenetic(phy))
#   phylomat[[i]] <- as.matrix(phydist)
# }
# 
# mat.new=matrix(nrow=nrow(com),ncol=nrow(com))
# for (k in 1:nrow(com)){
#   for (l in 1:nrow(com)){
#     
#     cm=numeric()
#     for(h in 1:nrow(phylo_res)){
#       cm=c(cm,phylomat[[h]][k,l])
#     }
#     mat.new[k,l]=mean(cm)
#   }
#   
# }
# 
# mean_phylo_mat<-data.frame(mat.new)
# row.names(mean_phylo_mat)<-row.names(phydist)
# names(mean_phylo_mat)<-names(phydist)
# 
# #write.csv(mean_phylo_mat, file ="mean_phylomat.csv")
# 
# #Sort mean_phylo_mat alphabetically by species name
# mean_phylo_mat2 <- (mean_phylo_mat[order(rownames(mean_phylo_mat)),order(colnames(mean_phylo_mat))])
# write.csv(mean_phylo_mat2, file="mean_phylomat2.csv")
# 
# 
# #create trait matrix
# threetraits<-traits[1:3]
# trait_mat_equal<-as.matrix(gowdis(threetraits))
# write.matrix(trait_mat_equal, file="mammal_trait_matrix.csv", sep=",")






############### FUNCTIONAL PHYLOGENETIC DISTANCE 

# focalSpecies <- data.frame(sp = c('Pecari_tajacu', 'Cervus_elaphus', 'Cuniculus_paca'))

# all combinations of species
focalSpecies <- c('Pecari_tajacu', 'Mazama_americana', 'Cuniculus_paca')

combos <- data.frame(Species1 = rep(focalSpecies, each = length(focalSpecies)),
                     Species2 = rep(focalSpecies, times = length(focalSpecies)))

# remove pairwise duplicates
# combos <- data.frame(t(apply(combos,1,sort)))
# combos <- combos[!duplicated(combos),]

# FDist <- FDist/max(FDist)
# PDist <- PDist/max(PDist)

# pull FDist and PDist from matrices to put into dataframe
combos$FPDist <- NA
for (i in 1:nrow(combos)){
  
  a <- 0.5 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
  p <- 2 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
  
  phylogeneticDistance <- PDist[combos[i,1], combos[i,2]]
  functionalDistance <- FDist[combos[i,1], combos[i,2]]
  
  combos$FPDist[i] <- ((a*(phylogeneticDistance^p)) + ((1-a)*(functionalDistance^p)))^(1/p)
  
}

# remove comparison between the same species 
combos <- combos[combos$Species1 != combos$Species2, ]
combos


# distance matrix with just focal species
FPDistMatrix <- matrix(ncol = length(focalSpecies), nrow = length(focalSpecies))

rownames(FPDistMatrix) <- focalSpecies
colnames(FPDistMatrix) <- focalSpecies

for (i in 1:nrow(FPDistMatrix)) {
  for (j in 1:ncol(FPDistMatrix)){
    
    phylogeneticDistance <- PDist[rownames(FPDistMatrix)[i], rownames(FPDistMatrix)[j]]
    functionalDistance <- FDist[rownames(FPDistMatrix)[i], rownames(FPDistMatrix)[j]]
    
    FPDistMatrix[i,j] <- ((a*(phylogeneticDistance^p)) + ((1-a)*(functionalDistance^p)))^(1/p)
    
  }
}


# distance matrix with all species
# all unique species 
allSpecies <- unique(c(colnames(FDist), colnames(PDist)))
FPDistMatrix <- matrix(NA, ncol = length(allSpecies), nrow = length(allSpecies))

rownames(FPDistMatrix) <- allSpecies
colnames(FPDistMatrix) <- allSpecies

# goal: make a matrix with FPDist for each species

for (i in 1:ncol(FPDistMatrix)){
  
  a <- 0.5 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
  p <- 2 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
  
  columnAnimal <- colnames(FPDistMatrix)[i]
  if (any(colnames(PDist) == columnAnimal)){ # make sure PDist has the animal in question
    if (any(colnames(FDist) == columnAnimal)) { # if both datasets have the animal in question
      
      
      
      for (j in 1:nrow(FPDistMatrix)){
        rowAnimal <- rownames(FPDistMatrix)[j]
        
        if (any(colnames(PDist) == rowAnimal)){ # make sure PDist has the animal in question
          if (any(colnames(FDist) == rowAnimal)) { # if both datasets have the animal in question
            
            phylogeneticDistance <- PDist[rowAnimal, columnAnimal]
            functionalDistance <- FDist[rowAnimal, columnAnimal]
            
            FPDistMatrix[rowAnimal, columnAnimal] <- ((a*(phylogeneticDistance^p)) + ((1-a)*(functionalDistance^p)))^(1/p)
          
            } else next
          } else next
        }  # closes j for loop
        
      
      
      
      } else next
    } else next
  }  # closes i for loop
  
  
  











for (i in 1:nrow(FPDistMatrix)) {
  rowAnimal <- rownames(FPDistMatrix)[i]
  if (i %% 5000){
    print(paste(i, "out of", nrow(FPDistMatrix, ":)")))
  }
  for (j in 1:ncol(FPDistMatrix)){
    columnAnimal <- colnames(FPDistMatrix)[j]
    a <- 0.5 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
    p <- 2 # based on https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12908
    
    if (any(colnames(PDist) == rowAnimal) & 
        any(colnames(PDist) == columnAnimal) &
        any(colnames(FDist) == rowAnimal) &
        any(colnames(FDist) == columnAnimal)){
      
      phylogeneticDistance <- PDist[rowAnimal, columnAnimal]
      functionalDistance <- FDist[rowAnimal, columnAnimal]
      
      FPDistMatrix[i,j] <- ((a*(phylogeneticDistance^p)) + ((1-a)*(functionalDistance^p)))^(1/p)
      
      
      
    } else next
    
  }
}
save(FPDistMatrix, file = 'allSpeciesFPDist.RData')













# make a list of biotic / FP distances for modeling

avgFPDist <- list()
for (i in 1:length(focalSpecies)) {
  distances <- FPDistMatrix[rownames(FPDistMatrix)[i],]
  avgDistanceFromOthers <- mean(distances[distances != 0], na.rm = TRUE)
  avgFPDist[[i]] <- rep(avgDistanceFromOthers, times = nrow(siteCovariate))
  # this will need to be edited when a list of species at each site is received
  # make new FPDist matrix for each site with only the species at that site
}

all(rownames(FPDistMatrix) == focalSpecies) # double check names are in correct order

names(avgFPDist) <- focalSpecies

save(avgFPDist, file = 'avgFPDist.RData')
save(FPDistMatrix, file = 'FPDistMatrix.RData')


