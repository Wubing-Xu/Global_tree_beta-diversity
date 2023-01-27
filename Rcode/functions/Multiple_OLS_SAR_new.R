
##A function to extract the statistical metrics of multiple OLS
OLScoef <- function(mydata){
  options(na.action="na.fail", use.fallback = TRUE)
  
  ## data preparation
  id<-!is.na(mydata[,1])
  mydata <- mydata[id, ]
  # full model formula
  name_pred <- colnames(mydata)[-1]
  n_pred <- length(name_pred)
  fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(name_pred, collapse= "+")))
  
  # the table for storing results
  result <- data.frame(matrix(NA, nrow = n_pred + 3, ncol=13))
  colnames(result) <- c("variable", "coef.full", "p.full", "coef.best", "coef.avg.full", "p.avg.full", 
                        "coef.avg.cod", "p.avg.cod", "coef_2.5", "coef_97.5", "weight", "R2", "uniqueR2")
  result[, 1] <- c(name_pred,"AIC","R2","weights")
  
  # fit the full model
  model_lm <- lm(fm, data=mydata)
  
  # get the coefficients of full model
  coefs <- summary(model_lm)$coef
  id <- match(rownames(coefs)[-1], result[,1])
  result[id, 2:3] <- coefs[-1, c(1, 4)]
  result[n_pred + 1, 2] <- AIC(model_lm)
  result[n_pred + 2, 2] <- round(summary(model_lm)$r.squared,3)
  
  # use multiple model inference
  mm_lm<- MuMIn::dredge(model_lm, rank="AIC")
  mm_coef <- as.data.frame(mm_lm)
  
  # the best model
  id <- match(colnames(mm_coef)[2:(n_pred + 1)], result[, 1])
  result[id, 4] <- as.numeric(mm_coef[1, 2:(n_pred + 1)]) # coefficients
  result[c(n_pred + 1, n_pred + 3), 4] <- as.numeric(mm_coef[1, c(n_pred + 4, n_pred + 6)]) # AIC and weight
  
  # the weight of full model
  id <- apply(mm_coef[, 2:(n_pred + 1)], 1, function(x) !any(is.na(x)))
  result[n_pred + 3, 2] <- mm_coef[id, n_pred + 6]
  
  # the averaged coef
  mm_avg <- model.avg(mm_lm)
  coef_ci <- confint(mm_avg)
  mm_avg <- summary(mm_avg)
  id <- match(rownames(mm_avg[[9]])[-1], result[,1])
  result[id, 5:6] <- mm_avg[[9]][-1, c(1, 5)]
  result[id, 7:8] <- mm_avg[[10]][-1, c(1, 5)]
  result[id, 9:10] <- coef_ci[-1, ]
  
  # the importance of each variable
  mm_imp <- MuMIn::sw(mm_avg)
  id <- match(names(mm_imp), result[, 1])
  result[id, 11] <- as.numeric(mm_imp)
  
  for(i in 1:n_pred){
    # calculate total R2 explained by each predictor variables
    subm_with <- as.formula(paste("~", paste(name_pred[i])))
    r2_sub_with <- summary(update(model_lm, subm_with))$r.squared
    result[i, 12] <- r2_sub_with
    
    # calculate unique R2 explained by each predictor variables
    subm_without <- as.formula(paste("~", paste(name_pred[-i], collapse= "+")))
    r2_sub_without <- summary(update(model_lm, subm_without))$r.squared
    result[i, 13] <- summary(model_lm)$r.squared - r2_sub_without
  }
  
  return(list("coef.table" = result, "fm_formula" = fm, "mm.lm" = mm_lm, "mm.avg" = mm_avg))
}


## A function to extract several statistical metrics using multiple SAR error model
MSARcoef <- function(coords, mydata, sardist = 150, longlat=TRUE){
  options(na.action="na.fail", use.fallback = TRUE)
  
  ## data preparation
  id<-!is.na(mydata[,1])
  mydata <- mydata[id, ]
  coords <- coords[id, ]
  # full model formula
  name_pred <- colnames(mydata)[-1]
  n_pred <- length(name_pred)
  fm <<- as.formula(paste(colnames(mydata)[1], "~", paste(name_pred, collapse= "+")))
  
  # the table for storing results
  result <- data.frame(matrix(NA, nrow = n_pred + 3, ncol=15))
  colnames(result) <- c("variable", "coef.full", "p.full", "coef.best", "coef.avg.full", "p.avg.full", 
                        "coef.avg.cod", "p.avg.cod", "coef_2.5", "coef_97.5", "weight", "R2", "uniqueR2", "ps.cor", "cor.p")
  result[, 1] <- c(name_pred,"AIC","R2","weights")
  
  # fit the full model
  nn <- dnearneigh(as.matrix(coords), d1=0, d2=sardist, longlat=longlat)  #or 200km
  w <- nb2listw(nn, zero.policy=TRUE)
  model_sar <- spatialreg::errorsarlm(fm, data=mydata, listw=w, zero.policy=T)
  
  # get the coefficients of full model
  coefs <- summary(model_sar)$Coef
  id <- match(rownames(coefs)[-1], result[,1])
  result[id, 2:3] <- coefs[-1, c(1, 4)]
  result[n_pred + 1, 2] <- AIC(model_sar)
  result[n_pred + 2, 2] <- as.numeric(cor(mydata[,1], predict(model_sar))^2) 
  
  # use multiple model inference
  mm_sar<- MuMIn::dredge(model_sar, rank="AIC")
  mm_coef <- as.data.frame(mm_sar)
  
  # the best model
  id <- match(colnames(mm_coef)[2:(n_pred + 1)], result[, 1])
  result[id, 4] <- as.numeric(mm_coef[1, 2:(n_pred + 1)]) # coefficients
  result[c(n_pred + 1, n_pred + 3), 4] <- as.numeric(mm_coef[1, c(n_pred + 4, n_pred + 6)]) # AIC and weight
  
  # the weight of full model
  id <- apply(mm_coef[, 2:(n_pred + 1)], 1, function(x) !any(is.na(x)))
  result[n_pred + 3, 2] <- mm_coef[id, n_pred + 6]
  
  # the averaged coef
  mm_avg <- model.avg(mm_sar)
  coef_ci <- confint(mm_avg)
  mm_avg <- summary(mm_avg)
  id <- match(rownames(mm_avg[[9]])[-c(1:2)], result[,1])
  result[id, 5:6] <- mm_avg[[9]][-(1:2), c(1, 4)]
  result[id, 7:8] <- mm_avg[[10]][-(1:2), c(1, 4)]
  result[id, 9:10] <- coef_ci[-(1:2), ]
  
  # the importance of each variable
  mm_imp <-  MuMIn::sw(mm_avg)
  id <- match(names(mm_imp), result[, 1])
  result[id, 11] <- as.numeric(mm_imp)
  
  # R2 of full model considering no spatial components of prediction
  fitted.full <- predict(model_sar, pred.type="trend")
  R2full <- cor(mydata[, 1], fitted.full)^2
  result[n_pred + 2, 12] <- R2full
  
  for(i in 1:n_pred){
    # calculate total R2 explained by each predictor variables
    subm_with <- as.formula(paste("~", paste(name_pred[i])))
    fitted_sub_with <- predict(update(model_sar, subm_with), pred.type="trend")
    result[i, 12] <- cor(mydata[, 1], fitted_sub_with)^2
    
    # calculate unique R2 explained by each predictor variables
    subm_without <- as.formula(paste("~", paste(name_pred[-i], collapse= "+")))
    fitted_sub_without <- predict(update(model_sar, subm_without), pred.type="trend")
    result[i, 13] <- R2full - cor(mydata[, 1], fitted_sub_without)^2
  }
  
  # use modified t-test to calculate pearson correlation
  for(i in 1:(n_pred)){
    z <- modified.ttest(x = unlist(mydata[, 1]), y = unlist(mydata[, i+1]), coords=coords, nclass=NULL)
    result[i, 14] <- round(z$corr, 4)
    result[i, 15] <- round(z$p.value, 4)
  }
  return(list("coef.table" = result, "fm_formula" = fm, "mm.sar" = mm_sar, "mm.avg" = mm_avg))
}


## A function to calculate R2 of bivariate OLS betwen one predictor with multiple response variables
OLS.mY.r2 <- function(x, y){
  result <- data.frame(matrix(NA, nrow=ncol(y), ncol=2))
  rownames(result) <- colnames(y)
  colnames(result) <- c("R2", "P")
  for(i in 1:ncol(y)){
    id <- !is.na(y[, i])
    xi <- x[id]
    yi <- y[id, i]
    model <- summary(lm(yi ~ xi))
    result[i, 1] <- round(model$adj.r.squared, 3)
    result[i, 2] <- coef(model)[2, 4]
  }
  return(result)
}


# A function to calculate R2 of bivariate OLS  betwen multiple predictor with one response variables
OLS.mX.r2 <- function(y, env){
  id <- !is.na(y)
  y <- y[id]
  env <- env[id,]
  
  result <- data.frame(matrix(NA, nrow=ncol(env), ncol=2))
  rownames(result) <- colnames(env)
  colnames(result) < -c("R2", "P")
  
  for(i in 1:ncol(env)){
    model <- summary(lm(y~env[,i]))
    result[i,1] <- round(model$adj.r.squared, 3)
    result[i,2] <- coef(model)[2,4]
  }
  return(result)
}


# A function to calculate correlation between multiple predictor with one response variables
cor.YmX <-function(y, env, method="pearson"){
  id <- !is.na(y)
  y <-y[id]
  env <- env[id,]
  
  result <- data.frame(matrix(NA, nrow=ncol(env), ncol=2))
  rownames(result) <- colnames(env)
  colnames(result) <- c("r", "P")
  
  for(i in 1:ncol(env)){
    model <- cor.test(y,env[,i],method=method)
    result[i,1] <- round(model$estimate,3)
    result[i,2] <- model$p.value
  }
  return(result)
}
