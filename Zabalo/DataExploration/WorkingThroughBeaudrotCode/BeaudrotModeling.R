######### R Code for Dynamic Occupancy Modeling #########
######### Model Set 2: Abiotic & Biotic         #########
######### Beaudrot, Acevedo, Lessard et al.     #########


rm(list=ls())
library(unmarked)
library(plyr)
library(AICcmodavg)


######################################
# Define helper functions
######################################

isEmpty <- function(x) {
  return(length(x)==0)
}


CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}


######################################
# Load formatted data
######################################


#load("All_covs.RData")
siteCovariate = read.csv("siteCovs2018.csv")
siteCovariate$Station <- as.factor(siteCovariate$Station)
siteCovariate$Hunting <- as.factor(siteCovariate$Hunting)
siteCovariate$Habitat <- as.factor(siteCovariate$Habitat)
siteCovariate$Community <- siteCovariate$Community/1000
siteCovariate$River <- siteCovariate$River/1000


#load("Species7sites_Include.RData")
peccary <- read.csv("CollaredPeccary.csv")
peccary <- data.frame(peccary[,-1], row.names=peccary[,1]) #rownames = stations

deer <- read.csv("Deer.csv")
deer <- data.frame(deer[,-1], row.names=deer[,1]) #rownames = stations

paca <- read.csv("Paca.csv")
paca <- data.frame(paca[,-1], row.names=paca[,1]) #rownames = stations

occupancyList <- list(Pecari_tajacu = peccary,
                      Mazama_americana = deer,
                      Cuniculus_paca = paca)


#load("BIOTIC_Include.RData")
load('avgFPDist.RData')


Species_data <- occupancyList
BIOTIC <- avgFPDist

nms = names(Species_data)

# Empty objects for loop
results.all=list()
mods.all=list()
results.table.ma=list()
results.table.aic=list()
colext.transformed=list()



############################################################################
################ Begin loop to run all models and extract results ##########
############################################################################

for(k in 1:length(nms)){
  print(k)
  
  # Define species for analysis and site using index value for list of all species
  sp.name <- names(Species_data)[k]
  sp.name
  species <- Species_data[[k]]
  site <- substr(names(Species_data)[[k]],1,3)
  site
  
  # Define covariates
  
  #siteCovs
  site_covs <- paste(site, "covs", sep="_")
  covs <- All_covs[names(All_covs)==site_covs]
  Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
  ForestLossCT <- unlist(as.matrix(sapply(covs, "[", 2)))
  ForestGainCT <- unlist(as.matrix(sapply(covs, "[", 3)))
  Biotic <- BIOTIC[[k]]
  site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic)
  
  #yearlySiteCovs
  Tmin <- as.data.frame(sapply(covs, "[", 4))
  Tmax <- as.data.frame(sapply(covs, "[", 5))
  Tvar <- as.data.frame(sapply(covs, "[", 6))
  Tsd <- as.data.frame(sapply(covs, "[", 7))
  Tmean <- as.data.frame(sapply(covs, "[", 8))
  
  # Define number of primary periods
  to=dim(Tmin)[2]
  years=as.character(1:to)
  years=matrix(years,nrow(species),to,byrow=TRUE)
  
  # Create object with data formatted for unmarked
  
  umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin, Tmax=Tmax, Tvar=Tvar, Tsd=Tsd, Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[2])
  
  # Create list to hold model set for model selection
  
  mods=list()
  
  # Null ##################################################################################
  try((fm0=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm0")) {
    if(CondNum(fm0)<2000){
      if(CondNum(fm0)>0){mods=c(mods,fm0)}
    } 
  }
  
  # Tmin only ##################################################################################
  try((fm2=colext(psiformula=~1,
                  gammaformula=~Tmin,
                  epsilonformula=~Tmin,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2")) {
    if(CondNum(fm2)<2000){
      if(CondNum(fm2)>0){mods=c(mods,fm2)}
    } 
  }
  
  try((fm2.1=colext(psiformula=~1,
                    gammaformula=~Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.1")) {
    if(CondNum(fm2.1)<2000){
      if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
    } 
  }
  
  try((fm2.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.2")) {
    if(CondNum(fm2.2)<2000){
      if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
    } 
  }
  
  
  # Tmax ##################################################################################
  try((fm3=colext(psiformula=~1,
                  gammaformula=~Tmax,
                  epsilonformula=~Tmax,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3")) {
    if(CondNum(fm3)<2000){
      if(CondNum(fm3)>0){mods=c(mods,fm3)}
    } 
  }
  
  try((fm3.1=colext(psiformula=~1,
                    gammaformula=~Tmax,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.1")) {
    if(CondNum(fm3.1)<2000){
      if(CondNum(fm3.1)>0){mods=c(mods,fm3.1)}
    }
  }
  
  
  try((fm3.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmax,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3.2")) {
    if(CondNum(fm3.2)<2000){
      if(CondNum(fm3.2)>0){mods=c(mods,fm3.2)}
    } 
  }
  
  # Tmean + Tsd only ##################################################################################
  try((fm4=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd,
                  epsilonformula=~Tmean + Tsd,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm4")) {
    if(CondNum(fm4)<2000){
      if(CondNum(fm4)>0){mods=c(mods,fm4)}
    }
  }
  
  try((fm4.1=colext(psiformula=~1,
                    gammaformula=~Tmean + Tsd,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm4.1")) {
    if(CondNum(fm4.1)<2000){
      if(CondNum(fm4.1)>0){mods=c(mods,fm4.1)}
    }
  }
  
  try((fm4.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmean + Tsd,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm4.2")) {
    if(CondNum(fm4.2)<2000){
      if(CondNum(fm4.2)>0){mods=c(mods,fm4.2)}
    }
  }
  
  
  
  # Tvar #############################################################################
  
  try((fm5=colext(psiformula=~1,
                  gammaformula=~Tvar,
                  epsilonformula=~Tvar,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5")) {
    if(CondNum(fm5)<2000){
      if(CondNum(fm5)>0){mods=c(mods,fm5)}
    } 
  }
  
  try((fm5.1=colext(psiformula=~1,
                    gammaformula=~Tvar,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.1")) {
    if(CondNum(fm5.1)<2000){
      if(CondNum(fm5.1)>0){mods=c(mods,fm5.1)}
    } 
  }
  
  try((fm5.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tvar,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.2")) {
    if(CondNum(fm5.2)<2000){
      if(CondNum(fm5.2)>0){mods=c(mods,fm5.2)}
    } 
  }
  
  
  # Biotic only ##################################################################################
  try((fm6=colext(psiformula=~1,
                  gammaformula=~Biotic,
                  epsilonformula=~Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6")) {
    if(CondNum(fm6)<2000){
      if(CondNum(fm6)>0){mods=c(mods,fm6)}
    }
  }
  
  try((fm6.1=colext(psiformula=~1,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.1")) {
    if(CondNum(fm6.1)<2000){
      if(CondNum(fm6.1)>0){mods=c(mods,fm6.1)}
    }
  }
  
  try((fm6.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm6.2")) {
    if(CondNum(fm6.2)<2000){
      if(CondNum(fm6.2)>0){mods=c(mods,fm6.2)}
    }
  }
  
  
  
  
  
  # Biotic + Tmin ##########################################################################
  try((fm7=colext(psiformula=~1,
                  gammaformula=~Biotic + Tmin,
                  epsilonformula=~Biotic + Tmin,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7")) {
    if(CondNum(fm7)<2000){
      if(CondNum(fm7)>0){mods=c(mods,fm7)}
    } 
  }
  
  try((fm7.1=colext(psiformula=~1,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.1")) {
    if(CondNum(fm7.1)<2000){
      if(CondNum(fm7.1)>0){mods=c(mods,fm7.1)}
    } 
  }
  
  try((fm7.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm7.2")) {
    if(CondNum(fm7.2)<2000){
      if(CondNum(fm7.2)>0){mods=c(mods,fm7.2)}
    } 
  }
  
  
  #Biotic * Tmin ##########################################################################
  try((fm8=colext(psiformula=~1,
                  gammaformula=~Biotic * Tmin,
                  epsilonformula=~Biotic * Tmin,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8")) {
    if(CondNum(fm8)<2000){
      if(CondNum(fm8)>0){mods=c(mods,fm8)}
    } 
  }
  
  try((fm8.1=colext(psiformula=~1,
                    gammaformula=~Biotic * Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.1")) {
    if(CondNum(fm8.1)<2000){
      if(CondNum(fm8.1)>0){mods=c(mods,fm8.1)}
    } 
  }
  
  try((fm8.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm8.2")) {
    if(CondNum(fm8.2)<2000){
      if(CondNum(fm8.2)>0){mods=c(mods,fm8.2)}
    } 
  }
  
  
  # Tmax + Biotic ##########################################################################
  try((fm9=colext(psiformula=~1,
                  gammaformula=~Tmax + Biotic,
                  epsilonformula=~Tmax + Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9")) {
    if(CondNum(fm9)<2000){
      if(CondNum(fm9)>0){mods=c(mods,fm9)}
    } 
  }
  
  try((fm9.1=colext(psiformula=~1,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.1")){
    if(CondNum(fm9.1)<2000){
      if(CondNum(fm9.1)>0){mods=c(mods,fm9.1)}
    } 
  }
  
  try((fm9.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm9.2")) {
    if(CondNum(fm9.2)<2000){
      if(CondNum(fm9.2)>0){mods=c(mods,fm9.2)}
    } 
  }
  
  
  
  
  #Tmax * Biotic ##########################################################################
  try((fm10=colext(psiformula=~1,
                   gammaformula=~Tmax * Biotic,
                   epsilonformula=~Tmax * Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10")) {
    if(CondNum(fm10)<2000){
      if(CondNum(fm10)>0){mods=c(mods,fm10)}
    } 
  }
  
  try((fm10.1=colext(psiformula=~1,
                     gammaformula=~Tmax * Biotic,
                     epsilonformula=~1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.1")) {
    if(CondNum(fm10.1)<2000){
      if(CondNum(fm10.1)>0){mods=c(mods,fm10.1)}
    } 
  }
  
  try((fm10.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmax * Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10.2")) {
    if(CondNum(fm10.2)<2000){
      if(CondNum(fm10.2)>0){mods=c(mods,fm10.2)}
    } 
  }
  
  
  # Tmean + Tsd + Biotic ##################################################################################
  try((fm11=colext(psiformula=~1,
                   gammaformula=~Tmean + Tsd + Biotic,
                   epsilonformula=~Tmean + Tsd + Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm11")) {
    if(CondNum(fm11)<2000){
      if(CondNum(fm11)>0){mods=c(mods,fm11)}
    }
  }
  
  try((fm11.1=colext(psiformula=~1,
                     gammaformula=~Tmean + Tsd + Biotic,
                     epsilonformula=~1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm11.1")) {
    if(CondNum(fm11.1)<2000){
      if(CondNum(fm11.1)>0){mods=c(mods,fm11.1)}
    }
  }
  
  
  try((fm11.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tmean + Tsd + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm11.2")) {
    if(CondNum(fm11.2)<2000){
      if(CondNum(fm11.2)>0){mods=c(mods,fm11.2)}
    }
  }
  
  
  
  # Tvar * Biotic ##########################################################################
  try((fm12=colext(psiformula=~1,
                   gammaformula=~Tvar * Biotic,
                   epsilonformula=~Tvar * Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12")) {
    if(CondNum(fm12)<2000){
      if(CondNum(fm12)>0){mods=c(mods,fm12)}
    } 
  }
  
  try((fm12.1=colext(psiformula=~1,
                     gammaformula=~Tvar * Biotic,
                     epsilonformula=~1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.1")) {
    if(CondNum(fm12.1)<2000){
      if(CondNum(fm12.1)>0){mods=c(mods,fm12.1)}
    } 
  }
  
  try((fm12.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar * Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.2")) {
    if(CondNum(fm12.2)<2000){
      if(CondNum(fm12.2)>0){mods=c(mods,fm12.2)}
    } 
  }
  
  
  
  
  #Tvar + Biotic ##########################################################################
  try((fm13=colext(psiformula=~1,
                   gammaformula=~Tvar + Biotic,
                   epsilonformula=~Tvar + Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13")) {
    if(CondNum(fm13)<2000){
      if(CondNum(fm13)>0){mods=c(mods,fm13)}
    } 
  }
  
  try((fm13.1=colext(psiformula=~1,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.1")) {
    if(CondNum(fm13.1)<2000){
      if(CondNum(fm13.1)>0){mods=c(mods,fm13.1)}
    } 
  }
  
  try((fm13.2=colext(psiformula=~1,
                     gammaformula=~1,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm13.2")) {
    if(CondNum(fm13.2)<2000){
      if(CondNum(fm13.2)>0){mods=c(mods,fm13.2)}
    } 
  }
  
  
  
  ######################################
  # Run Model Selection
  ######################################
  
  # If zero models for a species converged:
  
  if(isEmpty(mods)==TRUE){
    results.all[[k]]=NA  
    results.table.ma[[k]]=NA
    names(results.table.ma)[k] <- nms[k]
    results.table.aic[[k]]=NA
    names(results.table.aic)[k] <- nms[k]
    colext.transformed[[k]]=NA
    names(colext.transformed)[k] <- nms[k]
  }else{
    models=fitList(fits=mods)
    (ms <- modSel(models))
    results.all[[k]]=ms  	
    
    
    # For species with models that converged, run model selection:
    
    models=fitList(fits=mods)
    
    (ms <- modSel(models))
    
    results.all[[k]]=ms
    mods.all[[k]]=mods
    
    
    ######################################
    # Extract Results
    ######################################
    
    toExport<-as(ms,"data.frame")
    
    null.aic=toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]
    
    #if null didn't converge
    if(isEmpty(null.aic)==TRUE){
      null <- NA	
    }else{
      null <- toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",]
    }
    
    
    if((null.aic==0) ||isEmpty(null.aic)==TRUE){
      results.table.ma[[k]] <- rbind(null)
      temp=data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])}else{
        
        results.table.ma[[k]] <- rbind(toExport[1,],null)
        results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
        names(results.table.ma)[k] <- nms[k]
        
        temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
        names(temp) <- c("formula","delta","AICwt")
        results.table.aic[[k]] <- rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
        results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
        names(results.table.aic)[k] <- nms[k]
      }
    
    test=seq(3,length(toExport)-10,by=2)
    tmp <- as.numeric(toExport[1,test])
    
    colext.transformed[[k]] <- exp(tmp)
    colext.transformed[[k]] <- cbind(nms[k], toExport[1,2], colext.transformed[[k]])
    names(colext.transformed)[k] <- nms[k]
  }
  
  # Remove all models and results
  
  rm(fm0,
     fm2,fm2.1,fm2.2,
     fm3,fm3.1,fm3.2,
     fm4,fm4.1,fm4.2,
     fm5,fm5.1,fm5.2,
     fm6,fm6.1,fm6.2,
     fm7,fm7.1,fm7.2,
     fm8,fm8.1,fm8.2,
     fm9,fm9.1,fm9.2,
     fm10,fm10.1,fm10.2,
     fm11,fm11.1,fm11.2,
     fm12,fm12.1,fm12.2,
     fm13,fm13.1,fm13.2,
     mods,ms, tmp, temp, toExport)
}

############################################################################
# End loop
############################################################################




######################################
# Write Results
######################################

# Coerce the lists of results into dataframes and write to files

results.table.ma.df <- ldply(results.table.ma, data.frame)
results.table.aic.df <- ldply(results.table.aic, data.frame)
colext.transformed.df <- ldply(colext.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma.csv")
write.csv(results.table.aic.df, file="results.table.aic.csv")
write.csv(colext.transformed.df, file="colext.transformed.csv")
