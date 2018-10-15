---
title: "Multiclass classification"
date: 2018-01-28T21:55:52+01:00
anchor: "vowel"
weight: 42
---

Multiclass classification of Deterding's vowel data set. In each of the model runs, we checked whether the variational EM converged to different lower bounds, but it did not. So for the final model we do not require multiple restarts.

```R
# Load data set
data("vowel.train", package = "ElemStatLearn")
data("vowel.test", package = "ElemStatLearn")
vowel.train$y <- as.factor(vowel.train$y); n.train <- nrow(vowel.train)
vowel.test$y <- as.factor(vowel.test$y); n.test <- nrow(vowel.test)
vowel.dat <- rbind(vowel.train, vowel.test)

# Canonical I-probit
mod.can <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                   control = list(maxit = 10, restarts = TRUE))
iplot_par_lb(mod.can)
iplot_par_error(mod.can)
iplot_par_error(mod.can, type = "test")
mod.can <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                   control = list(maxit = 1000))

# fBm-0.5 I-probit
mod.fbm <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                   kernel = "fbm", control = list(maxit = 10, restarts = TRUE))
iplot_par_lb(mod.fbm)
iplot_par_error(mod.fbm)
iplot_par_error(mod.fbm, type = "test")
mod.fbm <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                   kernel = "fbm", control = list(maxit = 1000))


# SE I-probit
mod.se <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                  kernel = "se", control = list(maxit = 10, restarts = TRUE))
iplot_par_lb(mod.se)
iplot_par_error(mod.se)
iplot_par_error(mod.se, type = "test")
mod.se <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
                  kernel = "se", control = list(maxit = 1000))

# Confusion matrix
conf_mat <- function(mod) {
  tibble(test = vowel.test$y, predicted = mod$test$y) %>%
    group_by(test, predicted) %>%
    tally() %>%
    mutate(prop = n / 42) %>%
    complete(predicted, fill = list(n = 0, prop = 0)) -> plot.df

  levels(plot.df$test) <- levels(plot.df$predicted) <- c(
    "heed", "hid", "head", "had", "hud", "hard", "hod", "hoard", "hood", "who'd", "heard"
  )

  ggplot(plot.df, aes(test, predicted)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = ifelse(n > 0, n, NA)), na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(plot.df$test))) +
    scale_x_discrete(position = "top") +
    scale_fill_continuous(low = "white", high = "slategray", limits = c(0, 42)) +
    labs(x = "Test data", y = "Predicted classes") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0))
}

conf_mat(mod.can)
conf_mat(mod.fbm)
conf_mat(mod.se)

# Gaussian process classification using kernlab
# mod.gpc1 <- gausspr(y ~ ., vowel.train, kernel = "vanilladot")  # DOES NOT WORK
mod.gpc2 <- gausspr(y ~ ., vowel.train, kernel = "rbfdot")
y.hat <- predict(mod.gpc2, vowel.test)
RMSE.train <- sum(fitted(mod.gpc2) != vowel.train$y) / length(vowel.train$y) * 100
RMSE.test <- sum(y.hat != vowel.test$y) / length(vowel.test$y) * 100
```