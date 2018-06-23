---
title: "Comparison of methods"
date: 2018-01-28T21:55:52+01:00
anchor: "compareiprior"
weight: 37
---

Comparing direct optimisation (quasi-Newton methods), EM algorithm and MCMC for estimation of I-prior models.

```R
# Generate data set. Note that this data set can also be generated using
# iprior::gen_smooth().
set.seed(123)
N <- 150
f <- function(x, truth = FALSE) {
  35 * dnorm(x, mean = 1, sd = 0.8) +
    65 * dnorm(x, mean = 4, sd = 1.5) +
    (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
    3 * dnorm(x, mean = 2.5, sd = 0.3)
}
x <- c(seq(0.2, 1.9, length = N * 5 / 8), seq(3.7, 4.6, length = N * 3 / 8))
x <- sample(x, size = N)
x <- x + rnorm(N, sd = 0.65)  # adding random fluctuation to the x
x <- sort(x)
y.err <- rt(N, df = 1)
y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(N, mean = 4.1))  # adding random terms to the y

# Direct optimisation
mod1 <- iprior(y ~ x, dat.fit, kernel = "fbm", train.samp = 1:N)

# EM algorithm
mod2 <- iprior(y ~ x, dat.fit, kernel = "fbm", method = "em", train.samp = 1:N,
               control = list(maxit = 1000))

# MCMC 
stan.iprior.dat <- list(y = y, H = iprior::kern_fbm(x), n = length(y),
                        ytrue = f(x))
stan.iprior.mod <- "
data {
  int<lower=0> n; // number of data
  vector[n] y; // data y
  matrix[n, n] H; // the kernel matrix
  vector[n] ytrue; // true values
}
parameters {
  real alpha;
  real<lower=0> psi;
  real<lower=0> lambda;
}
transformed parameters {
  cov_matrix[n] Vy;
  Vy = psi * (lambda * H) * (lambda * H);
  for (i in 1:n)
    Vy[i, i] = Vy[i, i] + 1 / psi;
}
model {
  matrix[n, n] L_cov;
  L_cov = cholesky_decompose(Vy);
  y ~ multi_normal_cholesky(rep_vector(alpha, n), L_cov);
  lambda ~ normal(0, 10);
  psi ~ normal(0, 10);
}
generated quantities {
  // generate from posterior of y
  vector[n] w;
  matrix[n,n] Vw;
  vector[n] ynew;
  vector[n] muynew;
  matrix[n,n] Vynew;
  real rmse;
  real logLik;
  Vw = inverse(Vy);
  w = psi * lambda * H * inverse(Vy) * (y - alpha);
  muynew = alpha + lambda * H * w;
  Vynew = square(lambda) * H * Vw * H;
  for (i in 1:n)
    Vynew[i, i] = Vynew[i, i] + 1 / psi;
  ynew = multi_normal_rng(muynew, Vynew);
  rmse = sqrt(mean(square(ynew - ytrue)));
  logLik = multi_normal_lpdf(y|rep_vector(alpha, n), Vy);
}
"
# Compile the Stan programme
m <- stan_model(model_code = stan.iprior.mod)
m@model_name <- "iprior.fbm"

# Fit stan model
fit.stan <- sampling(m, data = stan.iprior.dat,
                     pars = c("alpha", "lambda", "psi", "ynew", "rmse", "logLik"),
                     iter = 2000, chains = 8, thin = 1)
print(fit.stan, pars = c("alpha", "lambda", "psi", "rmse", "logLik"))
```