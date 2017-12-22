## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
require(PenalizedAD)
data(bcg_data)

## ------------------------------------------------------------------------
f <- MAES(bcg_data, 2, attr(bcg_data, "label"))

## ------------------------------------------------------------------------
plot(f[[8]],col=attr(bcg_data, "label"),
     main="Classification performance",ylab="Predicted classes", xlab = "Observations")

## ------------------------------------------------------------------------
data(bcg_data)
data(experiment_data)

#1 fit the background as the potential background model if no signal is to be found
E_bcg=MAES(bcg_data,2,attr(bcg_data, "label"))

#fit the signal

#Iterate between fitting the background on data and fitting the signal on data2
minbic <- Inf
sig <- PAD(E_best=E_bcg, data_exp=experiment_data, gamma_eigen = 1,
            gamma_mean = 4,minbic=minbic,data=bcg_data)

#The component means
sig[[1]]

