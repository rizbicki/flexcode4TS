Code to implement the FlexCodeTS conditional density estimator described in the paper "Flexible conditional density estimation for time series", by G. Grivol, R. Izbicki, A. A. Okuno and R. B. Stern, available at https://arxiv.org/pdf/2301.09671.pdf

The method estimates $f(y_{t+1}|y_1,\ldots,y_t,x_t)$, where $x_t$ are exogenous variables using FlexCode (Izbicki, Rafael, and Ann B. Lee. "Converting high-dimensional regression to high-dimensional conditional density estimation." (2017): 2800-2831.)

Installation:

```r
install_github("flexcode4TS",username="rizbicki")
```

Example of usage:

```r
# Data generator to illustrate the method:
model_1 <- function(n)
{
  X <- data.frame(arima.sim(list(order=c(1,0,0), ar=.5), n=n))
  y <- rep(0,nrow(X))
  for(ii in 2:length(y))
  {
    y[ii] <- 0.5*y[ii-1] + X[ii,1] + 0.2*X[ii-1,1] + rnorm(1,0,0.9)
  }
  plot(y,type="l")
  return(list(X=X, y=y))
}

data <- model_1(n=10000)
# Calculate the number of training data points by subtracting 100 from the number of rows in 'X'
nTr <- nrow(data$X) - 100
# Split the data into training and new datasets based on 'nTr'
Xtrain <- data$X[1:nTr,,drop=FALSE]
ytrain <- data$y[1:nTr]
X_new <- data$X[-c(1:nTr),,drop=FALSE]
y_new <- data$y[-c(1:(nTr))]

# Fit FlexCodeTS to the time series data using the training set
fit <- fit_flexcode_timeseries(X=Xtrain, y=ytrain,
                               lags_x=3, lags_y=3,
                               nTrain=round(0.7*length(ytrain)),
                               regressionFunction=regressionFunction.XGBoost,
                               nIMax=20, chooseSharpen=TRUE)
# Plot the loss as a function of the number of expansion coefficients
plot(fit$cde_fit$errors, xlab="Number of expansion coefficients",
     ylab="Loss", cex.lab=1.38, pch=16)
# Print variable importance or other details from the fitted model
print(fit$cde_fit)

# Plot predictions for new data along with a prediction band at 95% probability
plot(fit, X_new, y_new, predictionBandProb=0.95,
     y_train = ytrain[(length(ytrain)-25):length(ytrain)])

``
