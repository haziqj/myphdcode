---
title: "Nyström example"
date: 2018-01-28T21:55:52+01:00
anchor: "nystrom"
weight: 50
---

Notet that fitting full model with `n = 2000` takes roughly 15 minutes.

```R
# Generate data set
dat <- gen_smooth(n = 2000, xlim = c(-1, 5.5), seed = 1)
head(dat)
ggplot2::qplot(X, y, data = dat)

# Full model
(mod.full <- iprior(y ~ X, dat, kernel = "fbm"))

# Nystrom method
(mod.nys <- iprior(y ~ X, dat, kernel = "fbm", nystrom = 50))

# Comparison and plots
get_time(mod.full); get_size(mod.full, "MB"); get_prederror(mod.full)
get_time(mod.nys); get_size(mod.nys); get_prederror(mod.nys)
plot(mod.full) + ggtitle("Full model")
plot(mod.nys) + ggtitle("Nyström model")
```