---
title: "Preliminary"
date: 2018-01-28T21:48:57+01:00
anchor: "preliminary"
weight: 15
---

Install and load all required `R` packages.

```r
if (!require("pacman")) install.packages("pacman")  # not strictly necessary, 
                                                    # but elegant way of loading 
                                                    # all packages
pacman::p_load(knitr, devtools, R2MLwiN, lme4, lmerTest, jmcm, caret, tidyverse,
	           directlabels, reshape2, kableExtra, cowplot, ElemStatLearn, 
	           spatstat, stringr, caret, kernlab, ggmcmc, coda, mvtnorm, 
	           MCMCglmm, runjags, rstan, plotly, parallel, doSNOW, RPEnsemble,
	           foreach, glmnet, BAS, ggmap, ggrepel)

# I-prior packages, install from GitHub
devtools::install_github("haziqj/iprior")
devtools::install_github("haziqj/iprobit")
devtools::install_github("haziqj/ipriorBVS")

# ggplot2 general setting
ggplot2::theme_set(theme_bw())

# rstan settings
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```