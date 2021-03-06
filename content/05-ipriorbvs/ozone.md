---
title: "Ozone data set"
date: 2018-01-28T21:55:52+01:00
anchor: "ozone"
weight: 44
---

Ozone data set (Breiman & Friedman and Casella & Moreno). Data set obtained from `mlbench` package.

```R
# Load data and clean up and create squared and interaction terms
data(Ozone, package = "mlbench")
colnames(Ozone) <- c("Month", "DayMonth", "DayWeek", "Ozone", "PresVand",
                     "WindLAX", "HumLAX", "TempSand", "TempElMon", "ibhLAX",
                     "PresGrad", "ibtLAX", "VisLAX")
Ozone <- Ozone[complete.cases(Ozone), ]
Ozone <- as.data.frame(lapply(Ozone, as.numeric))
y <- Ozone$Ozone; y <- scale(y)
X <- Ozone[, -4]; X <- scale(X)
X.sq <- X ^ 2
colnames(X.sq) <- paste0(colnames(X), "2")
X.int.names <- X.int <- NULL
for (i in seq_len(ncol(X))) {
  for (j in seq_len(ncol(X))) if (j > i) {
    X.int <- cbind(X.int, X[, j] * X[, i])
    X.int.names <- c(X.int.names, paste0(colnames(X)[j], "_", colnames(X)[i]))
  }
}
colnames(X.int) <- X.int.names
X2 <- cbind(X, X.sq, X.int)
n <- nrow(Ozone)

# This experiment uses 50 training points and the rest is used for testing 
# (obtaining out-of-sample test prediction error rates)
n.sim <- 50
res1 <- res2 <- data.frame(matrix(NA, nrow = n.sim, ncol = 4))
runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
for (i in seq_len(n.sim)) {
  test.samp <- sample(1:n, size = 25)
  train.samp <- (1:n)[-test.samp]

  y.test <- y[test.samp]
  y.train <- y[train.samp]
  X.test <- X[test.samp, ]
  X.train <- X[train.samp, ]
  X2.test <- X2[test.samp, ]
  X2.train <- X2[train.samp, ]

  mod1 <- ipriorBVS(y.train, X.train, two.stage = TRUE)
  mod2 <- ipriorBVS(y.train, X2.train, two.stage = TRUE)

  hpm1 <- get_pmps(mod1)[1]
  hpm2 <- get_pmps(mod2)[1]
  pred1 <- get_predict(mod1, X.test, y.test)
  pred2 <- get_predict(mod2, X2.test, y.test)

  res1[i, ] <- c(names(hpm1), as.numeric(hpm1), get_R2(mod1), sqrt(pred1$mse))
  res2[i, ] <- c(names(hpm2), as.numeric(hpm2), get_R2(mod2), sqrt(pred2$mse))

  cat("SIMULATION ", i, "\n")
}
colnames(res1) <- colnames(res2) <- c("Model", "PostProb", "R2", "RMSE")
res1[, -1] <- as.data.frame(lapply(res1[, -1], as.numeric))
res2[, -1] <- as.data.frame(lapply(res2[, -1], as.numeric))

res1 %>%
  group_by(Model) %>%
  summarise(prob = mean(`Post. Prob`), R2 = mean(R2), RMSE = mean(RMSE),
            prop = n() / 50) %>%
  arrange(desc(prob)) -> res1
res2 %>%
  group_by(Model) %>%
  summarise(prob = mean(`Post. Prob`), R2 = mean(R2), RMSE = mean(RMSE),
            prop = n() / 50) %>%
  arrange(desc(prob)) -> res2

# HPM for mod1
colnames(X)[as.logical(as.numeric(strsplit(res1[[1]][[1]], "")[[1]]))]

# HPM for mod2
colnames(X2)[as.logical(as.numeric(strsplit(res2[[1]][[1]], "")[[1]]))]
```