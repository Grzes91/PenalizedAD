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

