---
title: "Comparison of methods"
date: 2018-01-28T21:55:52+01:00
anchor: "compareiprobit"
weight: 38
---

Comparing Laplace approximation, variational EM and Hamiltonian MC for binary I-probit model. MCMC and Laplace approximation takes a very long time to complete.

```R
# Generate data set
dat <- iprobit::gen_spiral(n = 300, seed = 123)  # generate binary toy example data set

# Variational EM
mod.vi <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm")

# Laplace's method
mod.lap <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm",
                   method = "laplace")

# MCMC
stan.iprobit.dat <- list(y = as.numeric(dat$y) - 1, H = iprior::kern_fbm(dat$X),
                         n = length(dat$y))
stan.iprobit.mod <- "
data {
  int<lower=0> n; // number of data
  int<lower=0,upper=1> y[n]; // data y
  matrix[n, n] H; // the kernel matrix
}
parameters {
  real alpha;
  real<lower=0> lambda;
  vector[n] w;
}
transformed parameters {
  vector[n] mu;
  vector<lower=0,upper=1>[n] pi;
  mu = alpha + lambda * H * w;
  pi = Phi(mu);
}
model {
  y ~ bernoulli(pi);
  w ~ normal(0, 1);
  lambda ~ normal(0, 10);
}
generated quantities {
  // generate from posterior of y
  vector[n] ypred;
  vector[n] yvec;
  real brierscore;
  real errorrate;
  real logLik;
  for (i in 1:n)
    yvec[i] = y[i];
  brierscore = mean(square(pi - yvec));
  for (i in 1:n)
    if (mu[i] > 0) {
      ypred[i] = 1;
    } else {
      ypred[i] = 0;
    }
  errorrate = mean(square(ypred - yvec)) * 100;
  logLik = bernoulli_lpmf(y|pi);
}
"
# Compile the Stan programme
m <- stan_model(model_code = stan.iprobit.mod)
m@model_name <- "iprobit.fbm"

# Fit stan model
fit.stan <- sampling(m, data = stan.iprobit.dat,
                     pars = c("alpha", "lambda", "brierscore", "errorrate", "w",
                              "logLik"),
                     iter = 200, chains = 8, thin = 1)
print(fit.stan, pars = c("alpha", "lambda", "brierscore", "errorrate", "logLik"))

postmean <- summary(stan2coda(fit.stan))$stat[, 1]
postsd <- summary(stan2coda(fit.stan))$stat[, 2]
b.alpha <- postmean[grep("alpha", names(postmean))]
b.lambda <- postmean[grep("lambda", names(postmean))]
b.alpha.se <- postsd[grep("alpha", names(postsd))]
b.lambda.se <- postsd[grep("lambda", names(postsd))]
b.w <- postmean[grep("w", names(postmean))]
b.w.se <- postmean[grep("w", names(postsd))]

# Slight hack to get the picture
mod.hmc <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm",
                   control = list(theta0 = log(b.lambda), int.only = TRUE))
mod.hmc$w <- b.w
mod.hmc$param.full[1, ] <- b.alpha
iplot_predict(mod.hmc)
```