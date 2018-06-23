---
title: "I-prior BVS simulations"
date: 2018-01-28T21:55:52+01:00
anchor: "simulation"
weight: 40
---

There are two parts to the simulation study. The first part is immediately below, which involves all four of the priors/method: I-prior, g-prior, independent prior and also Lasso. This takes a long time to run, but for brevity the number of simulations is cut down from 100 to 4, so shouldn't take very long to test.

```R
# The simulation function to run all four methods in two stages simultaneously
other_res <- function(form, dat, method = c("lasso", "gprior"), truth) {
  method <- match.arg(method, c("lasso", "gprior"))
  if (method == "lasso") {
    X <- as.matrix(dat[, -1])
    Y <- as.numeric(dat$y)
    mod <- cv.glmnet(X, Y)
    where.lambdamin <- which(mod$lambda == mod$lambda.min)
    gamma <- as.numeric(mod$glmnet.fit$beta[, where.lambdamin] > 0)
    brier <- NA
  }
  if (method == "gprior") {
    p <- ncol(dat) - 1
    # First stage --------------------------------------------------------------
    mod <- bas.lm(as.formula(form), dat, prior = "g-prior", modelprior = uniform(),
                  alpha = p)
    tmp <- summary(mod)
    pips <- tmp[(1:p) + 1, 1]
    x.keep <- (pips >= 0.5)

    # Second stage -------------------------------------------------------------
    mod <- bas.lm(as.formula(form), dat[, c(TRUE, x.keep)], prior = "g-prior",
                  modelprior = uniform(), alpha = p)
    tmp <- summary(mod)
    gamma <- pips <- x.keep
    gamma[] <- pips[] <- 0
    pips[x.keep] <- tmp[seq_len(sum(x.keep)) + 1, 1]
    brier <- mean((pips - truth) ^ 2)
    gamma[x.keep] <- tmp[seq_len(sum(x.keep)) + 1, 2]
  }
  gamma.diff <- gamma - truth
  no.false.inclusions <- sum(gamma.diff > 0)  # false inclusions
  no.false.exclusions <- sum(gamma.diff < 0)  # false exclusions
  no.false.choices <- sum(gamma.diff != 0)  # false choices
  c(false.inc = no.false.inclusions, false.exc = no.false.exclusions,
    false = no.false.choices, brier = brier)
}

runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE, modules = "lecuyer")
bvs_sim <- function(snr = 0.90, n.sim = 100, save.intermediate = TRUE) {
  bvs.mods <- c("iprior_sing", "flat_prior", "gprior", "lasso")
  res <- matrix(NA, nrow = n.sim, ncol = 4)
  colnames(res) <- c("false.inc", "false.exc", "false", "brier")
  res <- list(res, res, res, res)  # one for each model
  names(res) <- bvs.mods

  for (sim in seq_len(n.sim)) {
    cat("\n")
    cat("------ SIMULATION NUMBER", sim, paste0("(SNR: ", snr, ")"), "------\n")
    cat("\n")
    dat <- gen_bvs(n = 150, p = 100, snr = snr)

    # Bayesian variable selection models ---------------------------------------
    for (i in seq_along(bvs.mods)) {  #seq_along(res)
      cat(bvs.mods[i])
      if (i < 3) {
        mod <- ipriorBVS(y ~ ., dat, model = bvs.mods[i], n.samp = 10000,
                         n.burnin = 500, n.adapt = 125, two.stage = TRUE
                         , priors = "lambda ~ dunif(0, 100)"
                         )
        tmp <- unlist(predict(mod, dat$truth))[1:4]
      } else {
        tmp <- other_res(y ~ ., dat$data, bvs.mods[i], dat$truth)
        cat("\n")
      }
      res[[i]][sim, ] <- tmp
      print(res[[i]][sim, ])
      cat("\n")
    }

    if (isTRUE(save.intermediate)) save(res, file = paste0("data/bvs-", snr))
  }

  push.message <- paste0(
    "BVS simulations (SNR: ", snr, ") COMPLETED."
  )
  pushoverr::pushover(message = push.message, user = userID, app = appToken)

  res
}

# Use the function above to do the simulations. Note that n.sim has been 
# reduced to 4 (originally 100) to reduce runtime if anyone wants to try it out
res.90 <- bvs_sim(snr = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.10, n.sim = 4)
res <- list("0.90" = res.90, "0.75" = res.75, "0.50" = res.50, "0.25" = res.25, "0.10" = res.10)

# Helper functions to do the analysis
count_false_choices <- function(x) {
  # Counts number of false choices
  x.0to2 <- mean(x <= 2)
  x.6 <- mean(x >= 6)
  x.3to5 <- 1 - x.0to2 - x.6
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

se_false_choices <- function(x) {
  # Gets the S.E. for number of false choices
  n <- length(x)
  x.0to2 <- (sd(x <= 2) * (n - 1) / n) / sqrt(n)
  x.6 <- (sd(x >= 6) * (n - 1) / n) / sqrt(n)
  x.3to5 <- (x > 2) & (x < 6)
  x.3to5 <- (sd(x.3to5) * (n - 1) / n) / sqrt(n)
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

bvs_res <- function(x = res, type = mean, models = "all") {
  # Get the results of BVS simulations. Choice of models are 
  # c("all", "iprior", "flat_prior", "gprior", "lasso")
  # type can be either "mean", "count_false_choices" or "se_false_choices"
  if (models == "all") {
    res <- lapply(x, function(y) t(sapply(y, function(z) apply(z, 2, type))))

    if (ncol(res[[1]]) > 3) {
      res <- lapply(res, function(y) {
        y <- y[, 7:9];
        colnames(y) <- c("0-2", "3-5", ">5");
        y
      })
    }
  } else {
    if (models == "iprior") tmp.res <- lapply(x, function(y) y[[1]])
    if (models == "flat_prior") tmp.res <- lapply(x, function(y) y[[2]])
    if (models == "gprior") tmp.res <- lapply(x, function(y) y[[3]])
    if (models == "lasso") tmp.res <- lapply(x, function(y) y[[4]])


    res <- t(sapply(tmp.res, function(x) apply(x, 2, type)))
    if (ncol(res) > 3) {
      res <- res[, 7:9]
      colnames(res) <- c("0-2", "3-5", ">5")
    }
  }
  res
}

# Run this to collate all results and plot
plot.iprior.df <- rbind(
  data.frame(res[[1]][[1]], SNR = "90%"),
  data.frame(res[[2]][[1]], SNR = "75%"),
  data.frame(res[[3]][[1]], SNR = "50%"),
  data.frame(res[[4]][[1]], SNR = "25%"),
  data.frame(res[[5]][[1]], SNR = "10%")
)

plot.flat.df <- rbind(
  data.frame(res[[1]][[2]], SNR = "90%"),
  data.frame(res[[2]][[2]], SNR = "75%"),
  data.frame(res[[3]][[2]], SNR = "50%"),
  data.frame(res[[4]][[2]], SNR = "25%"),
  data.frame(res[[5]][[2]], SNR = "10%")
)

plot.gprior.df <- rbind(
  data.frame(res[[1]][[3]], SNR = "90%"),
  data.frame(res[[2]][[3]], SNR = "75%"),
  data.frame(res[[3]][[3]], SNR = "50%"),
  data.frame(res[[4]][[3]], SNR = "25%"),
  data.frame(res[[5]][[3]], SNR = "10%")
)

plot.lasso.df <- rbind(
  data.frame(res[[1]][[4]], SNR = "90%"),
  data.frame(res[[2]][[4]], SNR = "75%"),
  data.frame(res[[3]][[4]], SNR = "50%"),
  data.frame(res[[4]][[4]], SNR = "25%"),
  data.frame(res[[5]][[4]], SNR = "10%")
)

ggplot(plot.iprior.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("I-prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p1

ggplot(plot.flat.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("Independent prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p2

ggplot(plot.gprior.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("g-prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p3

ggplot(plot.lasso.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("Lasso") +
  labs(x = "Number of false choices", y = "Proportion") -> p4

cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
```
