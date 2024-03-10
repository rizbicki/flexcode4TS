#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Create design matrix for time series analysis
#'
#' This function creates a design matrix from a given set of covariates (X) and a response variable (y). The design matrix includes lagged values of X and y.
#'
#' @param X a matrix or data frame of covariates
#' @param y a vector of response variable
#' @param lags_x number of lags for the covariates, defaults to NULL
#' @param lags_y number of lags for the response variable, defaults to the same value as lags_x
#'
#' @return X_new a matrix of the design matrix
create_design_matrix <- function(X=NULL,y,lags_x=NULL,
                                 lags_y=lags_x) {
  if(is.null(X)) # no covariates
  {
    if(lags_y==0)
      stop("need to use at least one lag for y if there are no covariates X")
    X_new <- matrix(Hmisc::Lag(y,shift=1))
    colnames(X_new)[ncol(X_new)] <-paste0("Y","_lag_",1)
    if(lags_y>1)
    {
      for(ii in 2:lags_y)
      {
        y_lag <- Hmisc::Lag(y,shift=ii)
        X_new <- cbind(X_new,y_lag)
        colnames(X_new)[ncol(X_new)] <-paste0("Y","_lag_",ii)
      }
      colnames(X_new)[1]<-paste0("Y","_lag_",1)
    }
  } else {
    if(!is.matrix(X))
      X <- as.matrix(X)
    X_new <- X
    colnames(X_new) <- paste0("X",1:ncol(X),"_lag_",0)
    if(lags_x>0)
    {
      for(ii in 1:lags_x)
      {
        X_lag <- apply(X,2,Hmisc::Lag,shift=ii)
        colnames(X_lag) <- paste0("X",1:ncol(X),"_lag_",ii)
        X_new <- cbind(X_new,X_lag)
      }
    }
    if(lags_y>0)
    {
      for(ii in 1:lags_y)
      {
        y_lag <- Hmisc::Lag(y,shift=ii)
        #if(is.vector(y_lag))
        #  y_lag <- t(y_lag)
        X_new <- cbind(X_new,y_lag)
        colnames(X_new)[ncol(X_new)] <-paste0("Y","_lag_",ii)
      }
    }
  }
  return(X_new)
}

#' Fits a FlexCoDE Model for Time Series Data
#'
#'
#' This function fits a FlexCoDE (Flexible Conditional Density Estimation) model for time series data. It is designed to handle input data \(X\) and response \(y\), incorporating lags for both \(X\) and \(y\). The function aims to train the model with a specified subset of observations and returns an object containing the model along with training data and specified lags.
#'
#' @param X A matrix or data frame containing the covariates for each observation. It represents the independent variables of the time series.
#' @param y A vector containing the response for each observation. It represents the dependent variable of the time series.
#' @param lags_x An integer specifying the number of lagged observations of \(X\) to include in the model. Default is NULL, indicating no lags unless specified.
#' @param lags_y An integer specifying the number of lagged observations of \(y\) to include in the model. If not specified, it defaults to the same number as `lags_x`.
#' @param nTrain An integer indicating the number of observations to be used for training. The default value is calculated as `round(0.8 * length(y))`, which uses 80% of the data for training.
#' @param ... Additional arguments that are passed directly to the underlying `fitFlexCoDE` function.
#
#' @return A list containing elements `lags_x`, `lags_y`, `X_train`, `y_train`, and
#'   `cde_fit`. The `cde_fit` element is the fitted model object from `fitFlexCoDE`.
#'   The return object has class `"fit_flexcode_timeseries"`.
#'
#'
#' @examples
#' # Generate sample time series data
#'
#' data <- generate_sample_data(n=10000)
#'
#' nTr <- nrow(data$X) - 100
#'
#' Xtrain <- data$X[1:nTr,, drop=FALSE]
#'
#' ytrain <- data$y[1:nTr]
#'
#' X_new <- data$X[-(1:nTr),, drop=FALSE]
#'
#' y_new <- data$y[-(1:nTr)]
#'
#' # Fit the FlexCoDE model
#'
#' fit <- fit_flexcode_timeseries(X=Xtrain, y=ytrain, lags_x=3, lags_y=3,
#'                                 nTrain=round(0.7 * length(ytrain)), regressionFunction=FlexCoDE::regressionFunction.XGBoost, nIMax=20,
#'                                 chooseSharpen=TRUE)
#'
#'
#' # Plot the model's errors
#'
#' plot(fit$cde_fit$errors, xlab="Number of expansion coefficients",
#'      ylab="Loss", cex.lab=1.38, pch=16)
#'
#' # Print model summary and plot variable importance
#'
#' print(fit$cde_fit)
#'
#' plot(fit, X_new, y_new, predictionBandProb=0.95,
#'      y_train=ytrain[(length(ytrain) - 25):length(ytrain)])
#'
#' @import FlexCoDE
#' @import ggplot2
#' @export
fit_flexcode_timeseries <-function(X=NULL,y,lags_x=NULL,
                                   lags_y=lags_x,nTrain=round(0.8*length(y)),
                                   ...)
{
  y <- as.vector(y)
  if(is.null(lags_x)&is.null(lags_y))
    stop("Provide at least one lag value for either x or y")

  return_value <- list(lags_x=lags_x,
                       lags_y=lags_y,
                       X_train=X,
                       y_train=y)

  X_design <- create_design_matrix(X=X,y=y,
                                   lags_x=lags_x,
                                   lags_y=lags_y)


  which_complete <- complete.cases(X_design)
  X_design <- X_design[which_complete,]
  y <- y[which_complete]

  if (is.null(ncol(X_design))){
    X_design = matrix(X_design,ncol=1)
  }

  #  random_index <- sample(1:length(y),
  #                         size = length(y),
  #                        replace = TRUE)
  random_index <- 1:length(y)
  fit <- fitFlexCoDE(xTrain = X_design[random_index[1:nTrain],],
                     zTrain=y[random_index[1:nTrain]],
                     xValidation = X_design[random_index[-c(1:nTrain)],],
                     zValidation=y[random_index[-c(1:nTrain)]],
                     ...)
  return_value$cde_fit=fit
  class(return_value) <- "fit_flexcode_timeseries"
  return(return_value)
}

#' @title predict_fit_flexcode_timeseries
#'
#' @description
#' predict_fit_flexcode_timeseries takes an object fit, a matrix X_new, a vector y_new, and a logical predictionBandProb as input. It returns the fitted model and the estimated density for the new observations.
#'
#' @param fit An object of class fit_flexcode_timeseries.
#' @param X_new A matrix containing the covariates for each new instance since the last observation used for training.
#' @param y_new A vector containing the response for each new instance since the last observation used for training. It has one last element than the number of rows of X_new. y_new must be NA if nrow(X_new)=1.
#' @param predictionBandProb A logical indicating whether to return the prediction band probability.
#'
#' @return A list with the predicted values and the estimated density.
#'
#' @examples
#' fit <- fit_flexcode_timeseries(X_train, y_train)
#' predict_fit_flexcode_timeseries(fit, X_new, y_new, predictionBandProb = TRUE)
#' @export
predict.fit_flexcode_timeseries <- function(fit,X_new=NULL,y_new=NA,
                                            predictionBandProb=FALSE)
{
  # X_new contains the covariates for each new instance since the last observation used for training
  # y_new contains the response for each new instance since the last observation used for training
  #   (it has one last element than the number of rows of X_new) y_new is NA if nrow(X_new)=1

  if(!is.na(y_new[length(y_new)]))
    y_new <- c(y_new,NA) # last y_new is the one that will be estimated

  if(is.null(X_new)) # no covariates
  {

  } else {

    if(nrow(X_new)==1)
    {
      if(!is.na(y_new))
        stop("y_new must be NA if X_new contains only one observation")
    }

    X_new_aug <- rbind(fit$X_train,X_new)

    y_new_aug <- c(fit$y_train,y_new)

    X_design_new <- create_design_matrix(X=X_new_aug,
                                         y=y_new_aug,
                                         lags_x=fit$lags_x,
                                         lags_y=fit$lags_y)
  }

  pred <- predict(fit$cde_fit,
                  X_design_new[nrow(X_design_new),,drop=FALSE],
                  predictionBandProb=predictionBandProb)
  pred_norm <- pred$CDE[1,]/sum(pred$CDE[1,])
  which_nn <- FNN::get.knnx(cumsum(pred_norm),
                            c((1-predictionBandProb)/2,
                              predictionBandProb+(1-predictionBandProb)/2),
                            k=1)$nn.index
  pred$pred_interval <- pred$z[which_nn]
  return(pred)
}

#' @title plot_fit_flexcode_timeseries
#'
#' @description
#' plot_fit_flexcode_timeseries takes an object fit, a matrix X_new, a vector y_new, a logical predictionBandProb, and a vector y_train as input and returns a plot of the estimated density for the new observations.
#'
#' @param fit An object of class fit_flexcode_timeseries.
#' @param X_new A matrix containing the covariates for each new instance since the last observation used for training.
#' @param y_new A vector containing the response for each new instance since the last observation used for training. It has one last element than the number of rows of X_new. y_new must be NA if nrow(X_new)=1.
#' @param predictionBandProb A logical indicating whether to return the prediction band probability.
#' @param y_train A vector containing the response for the training observations. If y_train is not provided, only the prediction region for the new observations will be plotted.
#'
#' @return A plot of the estimated density for the new observations.
#'
#' @export
plot.fit_flexcode_timeseries <- function(fit,X_new,
                                         y_new,predictionBandProb=TRUE,
                                         y_train=NULL)
{
  pred_values <- predict_experiments(fit,X_new=X_new,
                                     y_new=y_new,
                                     predictionBandProb = predictionBandProb)

  data=data.frame(x=pred_values$z,
                  y=pred_values$CDE[1,],
                  dataPoint=rep(1,length(pred_values$z)),
                  vertical=y_new[1])
  for(i in 2:length(y_new))
  {
    dataB=data.frame(x=pred_values$z,
                     y=pred_values$CDE[i,],
                     dataPoint=rep(i,length(pred_values$z)),
                     vertical=y_new[i])
    data=rbind(data,dataB)
  }
  g=ggplot(data,ggplot2::aes(x=x,y=y))+
    geom_line(size=1,color=2)+xlab("Response")+
    ylab("Estimated Density")+
    geom_vline(ggplot2::aes(xintercept=vertical),size=1)+
    theme(axis.title=ggplot2::element_text(size=2,face="bold"))+
    facet_wrap(~ dataPoint)+
    theme_bw()
  print(g)

  eps=0.35
  lineWidthPred=1.5
  k=length(y_new)

  if(!is.null(y_train))
  {
    plot(x=1:length(c(y_train,y_new)),y=c(y_train,y_new),
         main="",ylab="Prediction Region",
         cex.main=1.4,cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",
         xlim=c(0.5,length(c(y_train,y_new))+0.5),pch=16,
         ylim=c(fit$cde_fit$zMin,fit$cde_fit$zMax),xlab="Time",bty="l")
    for(ii in 1:k)
    {
      whichLarger=pred_values$CDE[ii,]>pred_values$th[ii]
      runs=rle(whichLarger>0)
      nRuns=length(runs$values)

      cumulative=cumsum(runs$lengths)
      for(jj in 1:nRuns)
      {
        if(runs$values[jj]==TRUE)
        {
          if(jj==1)
          {
            lower=fit$cde_fit$zMin
            upper=pred_values$z[cumulative[jj]]
            lines(c(ii+length(y_train),ii+length(y_train)),
                  c(lower,upper),col=1,lwd=lineWidthPred)
            lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
                  c(lower,lower),col=1,lwd=lineWidthPred)
            lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
                  c(upper,upper),col=1,lwd=lineWidthPred)
            next;
          }
          #points(rep(ii,sum(whichLarger)),predicted$z[whichLarger],pch=18,cex=0.9,col=2)
          lower=pred_values$z[cumulative[jj-1]]
          upper=pred_values$z[cumulative[jj]]
          lines(c(ii+length(y_train),ii+length(y_train)),
                c(lower,upper),col=1,lwd=lineWidthPred)

          lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
                c(lower,lower),col=1,lwd=lineWidthPred)
          lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
                c(upper,upper),col=1,lwd=lineWidthPred)
        }
      }
    }

    points(x=1:length(c(y_train,y_new)),y=c(y_train,y_new),main="",
           ylab="Estimate",cex.main=1.4,cex.axis=1.4,cex.lab=1.4,
           cex=1.5,col=2,xaxt="n",xlim=c(0.5,k+0.5),pch=16,
           xlab="Time")



    plot(x=1:length(c(y_train,y_new)),y=c(y_train,y_new),
         main="",ylab="Prediction Region",
         cex.main=1.4,cex.axis=1.4,cex.lab=1.4,cex=1.5,col=1,xaxt="n",
         xlim=c(0.5,length(c(y_train,y_new))+0.5),pch=16,
         ylim=c(fit$cde_fit$zMin,fit$cde_fit$zMax),xlab="Time",bty="l")
    for(ii in 1:k)
    {
      lower=pred_values$lower[ii]
      upper=pred_values$upper[ii]
      lines(c(ii+length(y_train),ii+length(y_train)),
            c(lower,upper),col=1,lwd=lineWidthPred)
      lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
            c(lower,lower),col=1,lwd=lineWidthPred)
      lines(c(ii-eps+length(y_train),ii+eps+length(y_train)),
            c(upper,upper),col=1,lwd=lineWidthPred)

    }
  }

  points(x=1:length(c(y_train,y_new)),y=c(y_train,y_new),main="",
         ylab="Estimate",cex.main=1.4,cex.axis=1.4,cex.lab=1.4,
         cex=1.5,col=2,xaxt="n",xlim=c(0.5,k+0.5),pch=16,
         xlab="Time")

  return()
}

predict_experiments <- function(fit,X_new=NULL,y_new=NULL,
                                predictionBandProb=0.95)
{
  # aux predict function for experiments (cv like test)


  pred_values <- predict(fit,X_new[1,,drop=FALSE],
                         predictionBandProb=predictionBandProb)
  estimated_densities <- matrix(NA,nrow(X_new),length(pred_values$z))
  th <- rep(NA,nrow(X_new))
  lower <- rep(NA,nrow(X_new))
  upper <- rep(NA,nrow(X_new))
  estimated_densities[1,] <- pred_values$CDE
  th[1] <- pred_values$th
  lower[1] <- pred_values$pred_interval[1]
  upper[1] <- pred_values$pred_interval[2]
  for(ii in 2:nrow(X_new))
  {
    pred_values <- predict(fit,X_new[1:ii,,drop=FALSE],
                           y_new = y_new[1:(ii-1)],
                           predictionBandProb=predictionBandProb)
    estimated_densities[ii,] <- pred_values$CDE
    th[ii] <- pred_values$th
    lower[ii] <- pred_values$pred_interval[1]
    upper[ii] <- pred_values$pred_interval[2]
  }
  return(list(z=pred_values$z,
              CDE=estimated_densities,
              th=th,lower=lower,
              upper=upper))
}
