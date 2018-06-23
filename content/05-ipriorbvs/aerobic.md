---
title: "Aerobic data set"
date: 2018-01-28T21:55:52+01:00
anchor: "aerobic"
weight: 42
---

Analysis of the Aerobic data set.

```R
# Load data from ipriorBVS package
data(aerobic, package = "ipriorBVS")
colnames(aerobic)[-1] <- paste0("X", 1:6)

# Fit one-stage model
(mod1 <- ipriorBVS(Oxygen ~ ., aerobic))
plot_coef2(mod1)

# Fit second-stage model
(mod2 <- ipriorBVS(Oxygen ~ ., aerobic, two.stage = TRUE))
```