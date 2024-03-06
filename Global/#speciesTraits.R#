########## Making a spreadsheet with our species and their traits

SNArecords <- read.csv("Siona/Data/SNAIndependentRecordsFormatted.csv") # just independent records
ZABrecords <- read.csv("Zabalo/Data/ZABIndependentRecordsFormatted.csv") # just independent records
SGErecords <- read.csv("Sinangoe/Data/SGEIndependentRecordsFormatted.csv") # just independent records
SKPrecords <- read.csv("Siekopai/Data/SKPIndependentRecordsFormatted.csv") # just independent records
records <- rbind(SNArecords[,c("Species","Station")], 
                 SGErecords[,c("Species","Station")], 
                 ZABrecords[,c("Species","Station")], 
                 SKPrecords[,c("Species","Station")])
sum(records$Species == "Mazama sp.")
records <- sort(unique(records$Species))
records <- records[! records %in% c("N/D N/D", "NAN NAN", "NA NA")] # remove N/D N/D, NAN NAN, NA NA
records <- data.frame(Name  = gsub(" ", "_", records),
                      Genus = sub(" .*", "", records),
                      Species = sub(".* ", "", records))
head(records)


pantheria <- read.delim("./FPDistData/PanTHERIA_1-0_WR05_Aug2008.txt", header = TRUE, sep = "\t") # https://esapubs.org/archive/ecol/E090/184/metadata.htm
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
traits$Name <- rownames(traits)
head(traits)
head(records)

traits[traits$Genus == "Mazama",]
together <- merge(records, traits[,c("Name","Order","Family")], by = "Name", all.x = TRUE, all.y = FALSE)
head(together)

write.csv(together, file = "Global/Data/speciesAttributes.csv")




