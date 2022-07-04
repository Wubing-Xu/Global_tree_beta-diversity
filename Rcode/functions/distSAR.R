## A function to do SAR analysis using different distance for selecitng a best distance

library(spdep)
library(ncf)

distSAR <- function(coords, y, env, start=101, end=1001, increment=100, co.increment=200, longlat=TRUE){
  
  # a data frame storing results
  dist_vect <- seq(start, end, by=increment)
  result<-data.frame("d" = c(0, dist_vect), "AIC"=NA, "Imax"=NA, "I1"=NA, "I20"=NA, "R2"=NA)
  
  # calculate the magnitudes of spatial autocorrelations in the residuals from OLS models  
  model_lm <- lm(y~.,data=env)
  co_lm <- correlog(coords[,1], coords[,2], resid(model_lm), increment=co.increment, resamp=0, latlon=longlat)
  result[1,2] <- AIC(model_lm)
  moranI20 <- co_lm$correlation[1:20]
  result[1,3] <- moranI20[which(abs(moranI20)==max(abs(moranI20)))] # maximum absolute value of Moran'r I at the first 20 distance classes
  result[1,4] <- moranI20[1] #  the value of Moran'r I at the first distance class
  result[1,5] <- sum(abs(moranI20)) # sum of the absolute value of Moran'r I at the first 20 distance classes
  result[1,6] <- as.numeric(cor(y, predict(model_lm))^2) # pseudo-R2, squared pearson correlation of predicted and observed values
  
  # calculate the magnitudes of spatial autocorrelations in the residuals from SAR error models.
  # a range of distances are used to calculate weight matrix
  for(i in 1:length(dist_vect)){
    # perform the SAR model
    nn <- dnearneigh(as.matrix(coords), d1=0, d2=dist_vect[i], longlat=longlat)
    w <- nb2listw(nn, style="W", zero.policy=TRUE)
    model_sar <- errorsarlm(y~., data=env, listw=w, zero.policy=T)
    
    # calculate the magnitudes of spatial autocorrelation 
    co_sar <- correlog(coords[,1], coords[,2], resid(model_sar), increment=co.increment, resamp=0, latlon=longlat)
    result[i+1,2] <- -2*model_sar$LL + 2*model_sar$parameters # AIC
    moranI20 <- co_sar$correlation[1:20]
    result[i+1,3] <- moranI20[which(abs(moranI20)==max(abs(moranI20)))]
    result[i+1,4] <- moranI20[1]
    result[i+1,5] <- sum(abs(moranI20))
    result[i+1,6] <- as.numeric(cor(y,predict(model_sar))^2)
  }
  
  return(result)
}
