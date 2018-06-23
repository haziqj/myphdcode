---
title: "Functional regression"
date: 2018-01-28T21:55:52+01:00
anchor: "tecator"
weight: 43
---

Regression with functional covariates (Tecator data set). Data obtained from `caret` package. Gaussian process regression models fitted using `kernlab`

```R
# Load data set
data(tecator, package = "caret")
endpoints <- as.data.frame(endpoints)
colnames(endpoints) <- c("water", "fat", "protein")

# Plot of 10 random spectra predicting fat content
# (Inspired by package caret, v6.0-68. See ?caret::tecator for details)
set.seed(123)
inSubset <- sample(1:dim(endpoints)[1], 10)

absorpSubset <- absorp[inSubset,]
endpointSubset <- endpoints$fat[inSubset]

newOrder <- order(absorpSubset[,1])
absorpSubset <- absorpSubset[newOrder,]
endpointSubset <- endpointSubset[newOrder]

as.tibble(data.frame(id = factor(1:10), fat = endpointSubset,
                     wave = absorpSubset)) %>%
  melt(id.vars = c("id", "fat")) %>%
  as.tibble() %>%
  mutate(x = rep(1:100, each = 10)) %>%
  ggplot(aes(x, value, col = id)) +
  geom_line() +
  geom_dl(aes(label = fat), method = list("last.polygons", cex = 0.9)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100)) +
  coord_cartesian(xlim = c(0, 105)) +
  labs(x = "Wavelength index", y = "Absorbance Units")

# Prepare data for iprior package
# n = 215, use first 172 for training, and remaining 43 for testing
absorp.train <- -t(diff(t(absorp)))	# this takes first differences using diff()
absorp.test <- absorp.train[173:215,]
absorp.train <- absorp.train[1:172,]

# Other variables
fat.train <- endpoints$fat[1:172]
fat.test <- endpoints$fat[173:215]
water.train <- endpoints$water[1:172]
water.test <- endpoints$water[173:215]

# NOTE: Some of the following model fits have set.seed() to fix the initial
# starting values for the estimation procedure. With direct optimisation,
# multiple local modes are a problem, but with EM it is less of a problem. All
# of the models were analysed with EM and multiple restarts, but for the
# vignette/paper, Models 1-5 showcase the iprior() functionality.

# Model 1: Canonical RKHS (linear)
set.seed(7) 
(mod1 <- iprior(y = fat.train, absorp.train))

# Model 2: Polynomial RKHS (quadratic)
set.seed(6) 
mod2 <- iprior(y = fat.train, absorp.train, kernel = "poly2",
               est.offset = TRUE)

# Model 3: Polynomial RKHS (cubic)
mod3 <- iprior(y = fat.train, absorp.train, kernel = "poly3",
               est.offset = TRUE)

# Model 4: fBm RKHS (default Hurst = 0.5)
(mod4 <- iprior(y = fat.train, absorp.train, kernel = "fbm",
                method = "em", control = list(stop.crit = 1e-3)))

# Model 5: fBm RKHS (estimate Hurst)
set.seed(17)
(mod5 <- iprior(fat.train, absorp.train, kernel = "fbm", method = "mixed",
                est.hurst = TRUE, control = list(stop.crit = 1e-3)))

# Model 6: SE kernel
(mod6 <- iprior(fat.train, absorp.train, est.lengthscale = TRUE,
                kernel = "se", control = list(restarts = TRUE,
                                              par.maxit = 100)))

# GPR models
mod7 <- gausspr(x = absorp.train, y = fat.train, kernel = "vanilladot")
y.hat <- predict(mod7, absorp.test)
RMSE.train7 <- sqrt(mod7@error)
RMSE.test7 <- sqrt(mean((y.hat - fat.test) ^ 2))

mod8 <- gausspr(x = absorp.train, y = fat.train)
y.hat <- predict(mod8, absorp.test)
RMSE.train8 <- sqrt(mod8@error)
RMSE.test8 <- sqrt(mean((y.hat - fat.test) ^ 2))
len.scale8 <- sqrt(1 / (2 * mod8@kernelf@kpar$sigma))

# Function to make it easier to compare models
tecator_res <- function(mod) {
  RMSE.train <- as.numeric(get_prederror(mod))
  RMSE.test <- sqrt(predict(mod, list(absorp.test), fat.test)$test.error)
  c("Log-lik" = logLik(mod), "Training RMSE" = RMSE.train,
    "Test RMSE" = RMSE.test)
}

tecator_tab <- function() {
  tab <- t(sapply(list(mod1, mod2, mod3, mod4, mod5, mod6),
                  tecator_res))
  rownames(tab) <- c("Linear", "Quadratic", "Cubic", "fBm-0.5",
                     paste0("fBm-", round(get_hurst(mod5), 3)),
                     paste0("SE-", round(get_lengthscale(mod6), 3)))
  tab
}

tecator_tab()
```