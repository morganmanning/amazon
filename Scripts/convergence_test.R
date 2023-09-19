
library(unmarked)

convergence_test <- function(data, period){
  
  y <- data[,2:49]
  
  closed_period = period
  
  ncol=48/closed_period
  
  mat1 <- matrix(0,ncol=ncol,nrow=30)
  
  starts <- seq(1,48,by=closed_period)
  limits <- seq(closed_period,48,by=closed_period)
  
  for (i in 1:30)
  {
    for (j in 1:ncol){
      if(all(is.na(y[i,starts[j]:limits[j]])) == TRUE){
        mat1[i,j]=NA
      } else if (sum(y[i,starts[j]:limits[j]],na.rm=TRUE) >= 1){
        mat1[i,j]=1
      }else {
        
        mat1[i,j]=0
      }
      
    }
  }
  
  umf <- unmarkedFrameOccu(y = mat1)
  
  fm <- occu(formula = ~1 ~1, data=umf)
  summary(fm)
  
  if(SE(fm)[1]=="NaN"){paste("Nope! It didn't converged")}else{
    paste("Congratulations! This model converged beautifully")
  }
  
}

#### Now test it
data_pec <- read.csv("CollaredPeccary.csv",header=TRUE)
data_deer <- read.csv("Deer.csv",header=TRUE)
data_paca <- read.csv("Paca.csv",header=TRUE)

convergence_test(data_paca, 4) # 4, 5, 7, 8, 9, 10, 13, 14, 15, 16




