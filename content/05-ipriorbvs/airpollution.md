---
title: "Mortality and air pollution"
date: 2018-01-28T21:55:52+01:00
anchor: "airpollution"
weight: 43
---

Mortality and air pollution data set (McDonald & Schwing's ridge analysis).

```R
# Load data set from iprior package
data("pollution", package = "iprior")
pollution.stand <- as.data.frame(scale(pollution))

# Fit the I-prior BVS two stage model
mod <- ipriorBVS(Mortality ~ ., pollution.stand, stand.x = FALSE, stand.y = FALSE,
                 two.stage = TRUE)
coef(mod)
```