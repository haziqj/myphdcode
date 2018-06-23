---
title: "Sensitivity analysis"
date: 2018-01-28T21:55:52+01:00
anchor: "sensitivity"
weight: 40
---


The second part involves only I-priors. This experiment is to see how well I-prior behaves under different SNR. Again, number of simulations have been reduced to 4 (originally 100) so as to provide a proof of concept.

```R
# Load the simulation functions
bvs_sim <- function(snr = 0.90, gamma = snr, n.sim = 100,
                    save.intermediate = TRUE) {
  res1 <- res2 <- matrix(NA, nrow = n.sim, ncol = 4)
  colnames(res1) <- colnames(res2) <- c("false.inc", "false.exc", "false", "brier")

  for (sim in seq_len(n.sim)) {
    cat("\n")
    cat("------ SIMULATION NUMBER", sim, paste0("(SNR: ", snr, ")"), "------")
    cat("\n")
    dat <- gen_bvs(n = 150, p = 100, snr = snr)

    # BVS I-prior (first stage) ------------------------------------------------
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE, modules = "lecuyer")
    mod <- ipriorBVS(y ~ ., dat, n.samp = 10000, n.burnin = 500,
                     n.adapt = 125, gamma.prob = gamma,
                     priors = "lambda ~ dunif(0, 100)", two.stage = FALSE)
    tmp <- unlist(predict(mod, dat$truth))[1:4]
    res1[sim, ] <- tmp
    cat("FIRST STAGE\n")
    print(res1[sim, ])
    cat("\n")

    # BVS I-prior (second stage) -----------------------------------------------
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE, modules = "lecuyer")
    mod <- ipriorBVS(y ~ ., dat, n.samp = 10000, n.burnin = 500,
                     n.adapt = 125, gamma.prob = gamma * get_mpm(mod),
                     priors = "lambda ~ dunif(0, 100)", two.stage = FALSE)
    tmp <- unlist(predict(mod, dat$truth))[1:4]
    res2[sim, ] <- tmp
    cat("SECOND STAGE\n")
    print(res2[sim, ])
    cat("\n")

    if (isTRUE(save.intermediate)) {
      save(res1, file = paste0("data/bvs-stage1-", snr, "-gam", gamma))
      save(res2, file = paste0("data/bvs-stage2-", snr, "-gam", gamma))
    }
  }

  push.message <- paste0(
    "BVS simulations (SNR: ", snr, ", gamma:", gamma, ") COMPLETED."
  )
  pushoverr::pushover(message = push.message, user = userID, app = appToken)

  list(stage1 = res1, stage2 = res2, snr = snr, gamma = gamma)
}

# This code runs the simulations for SNR = 0.90 under all gamma choices
res.90 <- bvs_sim(snr = 0.90, gamma = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.90, gamma = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.90, gamma = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.90, gamma = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.90, gamma = 0.10, n.sim = 4)
res.snr90 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])

# This code runs the simulations for SNR = 0.75 under all gamma choices
res.90 <- bvs_sim(snr = 0.75, gamma = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.75, gamma = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.75, gamma = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.75, gamma = 0.10, n.sim = 4)
res.snr75 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])

# This code runs the simulations for SNR = 0.50 under all gamma choices
res.90 <- bvs_sim(snr = 0.50, gamma = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.50, gamma = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.50, gamma = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.50, gamma = 0.10, n.sim = 4)
res.snr50 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])

# This code runs the simulations for SNR = 0.25 under all gamma choices
res.90 <- bvs_sim(snr = 0.25, gamma = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.25, gamma = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.25, gamma = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.25, gamma = 0.10, n.sim = 4)
res.snr25 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])

# This code runs the simulations for SNR = 0.10 under all gamma choices
res.90 <- bvs_sim(snr = 0.10, gamma = 0.90, n.sim = 4)
res.75 <- bvs_sim(snr = 0.10, gamma = 0.75, n.sim = 4)
res.50 <- bvs_sim(snr = 0.10, gamma = 0.50, n.sim = 4)
res.25 <- bvs_sim(snr = 0.10, gamma = 0.25, n.sim = 4)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.10, n.sim = 4)
res.snr10 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])

# Function to collate and plot the results
create_bvs_plot_df <- function(res, stage = 1) {
  rbind(
    data.frame(res[[1]][[stage]], gamma = "0.90"),
    data.frame(res[[2]][[stage]], gamma = "0.75"),
    data.frame(res[[3]][[stage]], gamma = "0.50"),
    data.frame(res[[4]][[stage]], gamma = "0.25"),
    data.frame(res[[5]][[stage]], gamma = "0.10")
  )
}

plot_res2 <- function(stage = 2, type = "false") {

  tmp <- function(res) {
    res <- matrix(rapply(res, function(x) apply(x, 2, mean, na.rm = TRUE)),
                  ncol = 4, byrow = TRUE)
    if (stage == 1) res <- res[c(1, 3, 5, 7, 9), ]
    if (stage == 2) res <- res[c(2, 4, 6, 8, 10), ]
    data.frame(res, gamma = c(0.90, 0.75, 0.50, 0.25, 0.10))
  }
  df <- data.frame(rbind(
    cbind(tmp(res.snr90), SNR = "90%"),
    cbind(tmp(res.snr75), SNR = "75%"),
    cbind(tmp(res.snr50), SNR = "50%"),
    cbind(tmp(res.snr25), SNR = "25%"),
    cbind(tmp(res.snr10), SNR = "10%")
  ))
  df <- df[, c(2, 1, 3, 4, 5, 6)]
  colnames(df)[1:4] <- c("false.exc", "false.inc", "false", "brier")
  df$SNR <- factor(df$SNR, levels = c("90%", "75%", "50%", "25%", "10%"))
  df <- reshape2::melt(df, id = c("gamma", "SNR"))
  levels(df$variable)[1:2] <- c("False exclusion", "False inclusion")

  ggplot(subset(df, df$var == "False inclusion" | df$var == "False exclusion"),
         aes(x = gamma, group = SNR)) +
    geom_line(aes(y = value,  col = SNR)) +
    facet_grid(. ~ variable) +
    labs(y = "Count", x = expression(paste("Hyperprior setting for ", pi[j]))) +
    theme_bw()
}
plot_res2()
```