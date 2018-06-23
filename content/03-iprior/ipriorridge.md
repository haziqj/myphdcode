---
title: "Ridge in log-likelihood"
date: 2018-01-28T21:55:52+01:00
anchor: "ridgeiprior"
weight: 39
---

Plot showing ridge in the I-prior log-likelihood. Uses the `plotly` package.

```R
# Fit an fBm I-prior model using the toy data set
mod <- iprior(y ~ X, gen_smooth(), kernel = "fbm")

# Create grid for plot
no.points <- 50
x <- log(get_lambda(mod))
x <- seq(x - 20, x + 15, length = no.points)  # lambda
x <- sort(c(x, log(get_lambda(mod))))
y <- log(get_psi(mod))
y <- seq(y - 10, y + 15, length = no.points)  # psi
y <- sort(c(y, log(get_psi(mod))))
tab <- expand.grid(x = x, y = y)

# Calculate the log-likelihood for each point in the grid
z <- rep(NA, nrow(tab))
for (i in seq_along(z)) {
  z[i] <- logLik(mod, theta = as.numeric(tab[i, ]))
}
tab.loglik <- matrix(z, nrow = no.points + 1, ncol = no.points + 1)
zmin <- -5000

# 3D plot
plot_ly(x = x, y = y, z = tab.loglik, zauto = FALSE, zmin = zmin) %>%
  add_surface() %>%
  layout(
    title = "I-prior log-likelihood",
    scene = list(
      xaxis = list(title = "log(lambda)"),
      yaxis = list(title = "log(psi)"),
      zaxis = list(title = "Log-likelihood",  range = c(zmin, max(z)))
    ))
```