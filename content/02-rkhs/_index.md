---
title: "Chapter 2"
date: 2018-01-28T22:01:36+01:00
anchor: "chapter2"
weight: 35
---

Plots of the sample I-prior paths in Chapter 2.

```R
# ggplot2 theme for the plots
theme_kernel_path <- theme_classic() +
  theme(legend.position = "none")

# RKHS of constant functions
kernel_path_constant <- function(n = 1000, seed = 789) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- matrix(1, nrow = n, ncol = n) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    coord_cartesian(ylim = c(min(f.prior) - 10, max(f.prior) + 10)) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample constant I-prior paths") +
    theme_kernel_path +
    scale_y_continuous(breaks = 0) +
    scale_x_continuous(breaks = 0)
}

# Canonical RKHS
kernel_path_canonical <- function(n = 1000, seed = 789, cen = TRUE,
                                  intercept = FALSE) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m))
    f.prior[, i] <- (as.numeric(intercept) +
                       iprior::kern_canonical(x, centre = cen)) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample linear I-prior paths") +
    theme_kernel_path +
    scale_y_continuous(breaks = 0) +
    scale_x_continuous(breaks = 0)
}

# Fractional Brownian motion RKHS
kernel_path_fbm <- function(n = 1000, seed = 789, hurst = 0.5, cen = TRUE,
                            intercept = FALSE) {
  the.title <- paste0("Sample fBm I-prior paths (Hurst = ", hurst, ")")
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m))
    f.prior[, i] <- (as.numeric(intercept) +
                       iprior::kern_fbm(x, gamma = hurst, centre = cen)) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = the.title) +
    theme_kernel_path +
    scale_y_continuous(breaks = 0) +
    scale_x_continuous(breaks = 0)
}

# Squared exponential RKHS
kernel_path_se <- function(n = 1000, seed = 789, l = 0.5, cen = TRUE,
                           intercept = FALSE) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m))
    f.prior[, i] <- (as.numeric(intercept) +
                       iprior::kern_se(x, l = l, centre = cen)) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = bquote(Sample~SE~"I-prior"~paths~(italic(l)~"="~.(l)))) +
    theme_kernel_path +
    scale_y_continuous(breaks = 0) +
    scale_x_continuous(breaks = 0)
}

# Pearson RKHS
kernel_path_pearson <- function(n = 1000, seed = 789) {
  x <- factor(sample(LETTERS[1:5], size = n, replace = TRUE))
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- iprior::kern_pearson(x) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- plot.df[match(unique(plot.df$x), plot.df$x), ]
  plot.df <- melt(plot.df, id = "x")
  ggplot(plot.df, aes(x, value)) +
    geom_point(aes(col = variable), size = 3, shape = 21, fill = NA, stroke = 1.2) +
    geom_point(aes(fill = variable, col = variable), size = 3, shape = 21, alpha = 0.5) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample Pearson I-prior points") +
    theme_kernel_path +
    scale_y_continuous(breaks = 0)
}

# Polynomial kernel RKHS
kernel_path_poly <- function(n = 1000, seed = 789, c = 0, d = 2, cen = TRUE,
                             intercept = 0.5) {
  the.title <- paste0("Sample polynomial I-prior paths (degree = ", d, ")")
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m))
    f.prior[, i] <- (as.numeric(intercept) +
                       iprior::kern_poly(x, c = c, d = d, centre = cen)) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = the.title) +
    theme_kernel_path +
    scale_y_continuous(breaks = 0) +
    scale_x_continuous(breaks = 0)
}
```
