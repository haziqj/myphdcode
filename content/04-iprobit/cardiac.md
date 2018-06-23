---
title: "Binary classification"
date: 2018-01-28T21:55:52+01:00
anchor: "cardiac"
weight: 40
---

Binary classification cardiac arrhythmia. Requires the arrhythmia data set `Arrh194.RData`. First load the required functions to run the simulation programme.

```R
# Set number of cores
no.cores <- parallel::detectCores()  

# Function to combine mean and se
mean_and_se <- function(x, y) paste0(dec_plac(x, 2), " (", dec_plac(y, 2), ")")

# Function to calculate ranks
tab_rank <- function(x, y) {
  # This is based on a weighted average. More weight given if low classification
  # error in small sammple size
  tmp <- apply(x + y, 1, function(x) sum(n * x) / sum(n))
  rank(tmp)
}

# Function to tabulate mean and se
tab_res <- function(...) {
  this <- list(...)
  K <- length(this)
  tab.mean <- tab.se <- tab <- NULL

  for (k in 1:K) {
    if (any(is.na(this[[k]]))) {
      tab.mean <- rbind(tab.mean, NA)
      tab.se <- rbind(tab.se, NA)
      tab <- rbind(tab, NA)
    } else {
      tab.mean.tmp <- apply(this[[k]], 2, mean)
      tab.se.tmp <- apply(this[[k]], 2, sd) / sqrt(nrow(this[[k]]))
      tab.mean.and.se <- mean_and_se(tab.mean.tmp, tab.se.tmp)
      tab.mean <- rbind(tab.mean, tab.mean.tmp)
      tab.se <- rbind(tab.se, tab.se.tmp)
      tab <- rbind(tab, tab.mean.and.se)
    }
  }

  rownames(tab.mean) <- rownames(tab.se) <- rownames(tab) <- names(this)
  colnames(tab.mean) <- colnames(tab.se) <- colnames(tab) <- colnames(this[[1]])
  list(
    tab = as.data.frame(tab),
    tab.mean = as.data.frame(tab.mean),
    tab.se = as.data.frame(tab.se)
  )
}

# I-probit simulation
my_iprobit_sim <- function(nsim = 100, kernel = c("canonical", "fbm", "se")) {
  # kernel <- match.arg(kernel, c("canonical", "fbm", "se"))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)

  cl <- makeCluster(no.cores)
  registerDoSNOW(cl)
  res <- foreach(i = 1:nsim, .combine = rbind,
                 .packages = c("iprior", "iprobit"),
                 .export = c("n", "y", "X"),
                 .options.snow = list(progress = progress)) %dopar% {
      res.tmp <- rep(NA, length(n))
      for (j in 1:length(n)) {
        tmp <- iprobit(y, X, kernel = kernel, control = list(maxit = 5),
                       train.samp = sample(seq_along(y), size = n[j]))
        res.tmp[j] <- predict(tmp)$error
      }
      res.tmp
    }
  close(pb)
  stopCluster(cl)

  push.message <- paste0(
    experiment.name, ": ", kernel, " I-prior probit COMPLETED."
  )
  pushoverr::pushover(message = push.message, user = userID, app = appToken)

  colnames(res) <- paste0(c("n = "), n)
  res
}

# 451 observations with 194 continuous covariates
load("data/Arrh194.RData")
experiment.name <- "Cardiac data"
X.orig <- ArrhDataNew$x
y <- ArrhDataNew$y
y <- y - 1  # convert to 0 and 1
y <- as.factor(y)
levels(y) <- c("Normal", "Arrhythmia")
N <- length(y)
n <- c(50, 100, 200)  # subsamples
X <- scale(X.orig)  # standardise data
X[, 36] <- rep(0, N); X[, 181] <- rep(0, N)  # all zeroes
summary(y)

# Full I-probit model
set.seed(123)
(mod <- iprobit(y, X, kernel = "fbm", control = list(maxit = 100)))
iplot_lb(mod, lab.pos = "down") +
  scale_y_continuous(sec.axis = dup_axis(name = " ")) +
  coord_cartesian(xlim = c(1, 14))
iplot_error(mod) + coord_cartesian(xlim = c(1, 14))

# Simulation study
set.seed(456)
res.canonical <- my_iprobit_sim(nsim = 100, kernel = "canonical")
res.fbm       <- my_iprobit_sim(nsim = 100, kernel = "fbm")
res.se        <- my_iprobit_sim(nsim = 100, kernel = "se")

(tab <- tab_res("I-probit (linear)"  = res.canonical,
                "I-probit (fBm-0.5)" = res.fbm,
                "I-probit (SE-1.0)"  = res.se))

# Other results
knn.mean        <- c(40.64, 38.94, 35.76)
knn.se          <- c(0.33, 0.33, 0.36)
knn             <- mean_and_se(knn.mean, knn.se)
svm.mean        <- c(36.16, 35.64, 35.20)
svm.se          <- c(0.47, 0.39, 0.35)
svm             <- mean_and_se(svm.mean, svm.se)
svm.radial.mean <- c(48.39, 47.24, 46.85)
svm.radial.se   <- c(0.49, 0.46, 0.43)
svm.radial      <- mean_and_se(svm.radial.mean, svm.radial.se)
gpc.radial.mean <- c(37.28, 33.80, 29.31)
gpc.radial.se   <- c(0.42, 0.40, 0.35)
gpc.radial      <- mean_and_se(gpc.radial.mean, gpc.radial.se)
rf.mean         <- c(31.65, 26.72, 22.40)
rf.se           <- c(0.39, 0.29, 0.31)
rf              <- mean_and_se(rf.mean, rf.se)
nsc.mean        <- c(34.98, 33.00, 31.08) # nearest shrunken centroids
nsc.se          <- c(0.46, 0.440, 0.41)
nsc             <- mean_and_se(nsc.mean, nsc.se)
penlog.mean     <- c(34.92, 30.48, 26.12) # L1-penalised logistic regression
penlog.se       <- c(0.42, 0.34, 0.27)
penlog          <- mean_and_se(penlog.mean, penlog.se)

other.tab <- rbind(
  "k-nn"           = knn,
  "SVM"            = svm,
  "GP (radial)"    = gpc.radial,
  "Random forests" = rf,
  "NSC"            = nsc,
  "L-1 logistic"   = penlog
)
colnames(other.tab) <- colnames(tab$tab)

# Calculate ranks
tab.mean <- rbind(tab$tab.mean,
                  "k-nn"           = knn.mean,
                  "SVM"            = svm.mean,
                  "GP (radial)"    = gpc.radial.mean,
                  "Random forests" = rf.mean,
                  "NSC"            = nsc.mean,
                  "L-1 logistic"   = penlog.mean
)
tab.se <- rbind(tab$tab.se,
                "k-nn"           = knn.se,
                "SVM"            = svm.se,
                "GP (radial)"    = gpc.radial.se,
                "Random forests" = rf.se,
                "NSC"            = nsc.se,
                "L-1 logistic"   = penlog.se
)
tab.ranks <- tab_rank(tab.mean, tab.se)

# Tabulate results
cbind(rbind(tab$tab, other.tab), Rank = tab.ranks)
```