---
title: "Spatio-temporal modelling"
date: 2018-01-28T21:55:52+01:00
anchor: "btb"
weight: 43
---

Spatio-temporal modelling of bovine tuberculosis (BTB) in Cornwall. Requires `BTBppp.RData` and `BTB_spoligotype_data.txt`. Models obviously take a long time to fit, but so do the plotting functions.

```R
# Load data set
load("data/BTBppp.RData")
plot.df <- as.data.frame(pppdata$window)
W <- pppdata$window
simpW <- simplify.owin(W, 100)  # This reduces resolution of plot
plot.df <- as.data.frame(W)

# Get data points (location, spoligotypes and year)
btb <- read.table("data/BTB_spoligotype_data.txt", header = FALSE, sep = " ",
                  skip = 1)
colnames(btb) <- c("x", "y", "year", "sp")
btb$sp <- factor(btb$sp)

# Keep the four largest classes
btb %>%
  group_by(sp) %>%
  tally() %>%
  arrange(desc(n)) -> btb.n
(classes.to.keep <- as.character(btb.n$sp[1:4]))
levels(btb$sp) <- c("Sp9",    "Others", "Others", "Sp12" ,  "Sp15",
                    "Others", "Sp20",   "Others", "Others", "Others")
btb$sp <- factor(btb$sp, levels = c("Sp9", "Sp12", "Sp15", "Sp20", "Others"))

# A similar table to Diggle et al. (2005)
btb.summary <- btb %>%
  group_by(year, sp) %>%
  arrange(year) %>%
  tally() %>%
  complete(sp, year, fill = list(n = 0)) %>%
  reshape2::dcast(year ~ sp)

# Group the years into 4 categories: 1) < 1997; 2) 1997-1998; 3) 1999-2000; 4)
# 2001-02
btb_year_fn <- function(x) {
  res <- x
  res[x < 1997] <- "< 1997"
  res[x == 1997 | x == 1998] <- "1997-1998"
  res[x == 1999 | x == 2000] <- "1999-2000"
  res[x == 2001 | x == 2002] <- "2001-2002"
  res
}
as.tibble(btb) %>%
  mutate(period = factor(btb_year_fn(year))) -> btb

# Bar plot
as.tibble(btb.summary) %>%
  reshape2::melt(id.vars = "year") %>%
  ggplot(aes(year, value, col = variable)) +
  geom_segment(aes(xend = year, yend = 0), size = 12) +
  scale_x_continuous(breaks = seq(1989, 2002, by = 1),
                     minor_breaks = seq(1989, 2002, by = 1)) +
  labs(x = "Year", y = "Count") +
  scale_color_discrete(name = "Spoligotype") +
  theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 6)))

# Plot of Cornwall
ggplot() +
  geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA, col = "grey25") +
  geom_point(data = btb, aes(x, y, col = sp)) +
  labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)",
       col = "Spoligotype") +
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(labels = function(x) x / 1000) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.005), legend.justification = c(1, 0))

# Prepare to fit I-probit models
Spoligotype <- btb$sp
X <- scale(btb[, 1:2])
mu.x <- attr(X, "scaled:center")
sd.x <- attr(X, "scaled:scale")
year <- scale(btb$year, scale = FALSE)
mu.year <- attr(year, "scaled:center")
period <- btb$period

# Constant only model
mod0 <- iprobit(y = Spoligotype, X, control = list(int.only = TRUE, theta0 = -100))

# Spatial model
mod1 <- iprobit(y = Spoligotype, X, kernel = "fbm",
                control = list(maxit = 200, restarts = TRUE))

# Spatio-period model (period is categorical)
mod2 <- iprobit(y = Spoligotype, X, period, kernel = "fbm",
                interactions = "1:2", control = list(maxit = 200, restarts = TRUE))

# Spatio-temporal model (year is continuous)
mod3 <- iprobit(y = Spoligotype, X, year, kernel = "fbm",
                interactions = "1:2", control = list(maxit = 200, restarts = TRUE))

# Tabulate results
rest_fn <- function(mod) {
  res <- iprior::dec_plac(summary(mod)$tab[, c(1, 2)], 3)

  mod.name <- deparse(substitute(mod))
  if (mod.name == "mod0")
    res <- rbind(res[-6, ], NA, NA)
  if (mod.name == "mod1")
    res <- rbind(res, NA)

  pij <- fitted(mod)$p
  suppressWarnings(y <- iprior::get_y(mod$ipriorKernel))
  loglik <- sum(log(pij[cbind(seq_along(y), y)]))  # matrix indexing

  rbind(
    res,
    # loglik = c(iprior::dec_plac(logLik(mod), 2), NA),
    loglik = c(iprior::dec_plac(loglik, 2), NA),
    error  = c(iprior::dec_plac(get_error_rate(mod), 2), NA),
    brier  = c(iprior::dec_plac(get_brier_score(mod), 3), NA)
  )
}
tab.btb <- cbind(rest_fn(mod0), " " = "",
                 rest_fn(mod1), " " = "",
                 rest_fn(mod3), " " = "",
                 rest_fn(mod2))
rownames(tab.btb) <- c(
  paste0("Intercept (", levels(btb$sp), ")"),
  "Scale (spatial)", "Scale (temporal)",
  "Log-likelihood", "Error rate (%)", "Brier score"
)
tab.btb

# Preparation for plotting the predicted probability surfaces. First
# obtain points inside the polygon (takes a while)
rescalee <- function(x) {
  res <- x * rep(sd.x, each = nrow(x)) + rep(mu.x, each = nrow(x))
  colnames(res) <- c("x", "y")
  res
}
X.var <- 1:2
maxmin <- cbind(apply(X, 2, min), apply(X, 2, max))
xx <- list(NULL)
for (j in 1:2) {
  mm <- maxmin[X.var[j], ]
  xx[[j]] <- seq(from = mm[1] - 1, to = mm[2] + 1, length.out = 300)  # change this for resolution
}
mm <- maxmin[X.var, ]
x.df.full <- expand.grid(xx[[1]], xx[[2]])
tmp <- x.df.full * rep(sd.x, each = nrow(x.df.full)) +
  rep(mu.x, each = nrow(x.df.full))
isin <- inside.owin(tmp[, 1], tmp[, 2], W)
x.df <- x.df.full[isin, ]

# Calculate fitted probabilities (also takes a while!)
fill.col <- iprior::gg_col_hue(5)
N <- nrow(x.df)
a <- predict(mod1, list(x.df))  # spatial-model

period.df <- data.frame(
  factor(rep(levels(btb$period)[1], N), levels = levels(btb$period)),
  factor(rep(levels(btb$period)[2], N), levels = levels(btb$period)),
  factor(rep(levels(btb$period)[3], N), levels = levels(btb$period)),
  factor(rep(levels(btb$period)[4], N), levels = levels(btb$period))
)

a1 <- predict(mod2, list(x.df, period.df[, 1]))
a2 <- predict(mod2, list(x.df, period.df[, 2]))
a3 <- predict(mod2, list(x.df, period.df[, 3]))
a4 <- predict(mod2, list(x.df, period.df[, 4]))

b <- list(NULL)
unq.year <- sort(as.numeric(unique(year)))
for (i in seq_along(unique(year))) {
  b[[i]] <- predict(mod3, list(x.df, matrix(unq.year[i], nrow = N, ncol = 1)))$prob
  print(i)
}

# Function to plot the SPATIAL MODEL
plot_spatial_model <- function(m = 1, method = "bottom.pieces", points = FALSE,
                               contour.labels = FALSE) {
  current.label <- paste0("Spoligotype ", gsub("Sp", "", levels(btb$sp)[m]))
  contour.df <- cbind(rescalee(x.df), prob = a$prob[, m])
  ggplot(data = contour.df, aes(x, y)) +
    geom_raster(aes(fill = prob)) -> v

  v <- v + geom_contour(aes(z = prob), col = "grey20", size = 0.3) +
    geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA,
                 col = "grey25") +
    labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_continuous(name = NULL, low = "white", high = fill.col[m],
                          limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_bw() +
    annotate("text", x = min(plot.df$x), y = max(plot.df$y),
             label = current.label, vjust = 1, hjust = 0) +
    theme(panel.grid = element_blank())

  if (isTRUE(points)) {
    v <- v + geom_point(data = subset(btb, sp == levels(btb$sp)[m]), aes(x, y),
                        col = "grey30", shape = 21, fill = fill.col[m], size = 1)
  }

  if (isTRUE(contour.labels)) {
    v <- v + directlabels::geom_dl(aes(z = prob, label = ..level..),
                                   stat = "contour", colour = "grey20",
                                   method = list("far.from.others.borders",
                                                 "calc.boxes", "draw.rects",
                                                 cex = 0.65)) +
      guides(fill = FALSE)
  } else {
    v <- v + guides(fill = guide_colorbar(barwidth = 0.45, barheight = 16))
  }

  v
}

p1 <- plot_spatial_model(1, points = TRUE, contour.labels = TRUE)
p2 <- plot_spatial_model(2, points = TRUE, contour.labels = TRUE)
p3 <- plot_spatial_model(3, points = TRUE, contour.labels = TRUE)
p4 <- plot_spatial_model(4, points = TRUE, contour.labels = TRUE)
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2,
                   labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 10,
                   label_fontface = "plain")

# Function to plot for SPATIO-TEMPORAL MODEL
mean.year <- attr(year, "scaled:center")
plot_stemporal_model <- function(year = 1, points = TRUE,
                                 mod = c("period", "temporal")) {
  mod <- match.arg(mod, c("period", "temporal"))

  if (mod == "period") {
    current.period <- levels(btb$period)[year]
    current.label <- paste0("Year: ", current.period)
    alphaa <- 0.95
    if (year == 1)
      contour.df <- cbind(rescalee(x.df), prob = a1$prob)
    if (year == 2)
      contour.df <- cbind(rescalee(x.df), prob = a2$prob)
    if (year == 3)
      contour.df <- cbind(rescalee(x.df), prob = a3$prob)
    if (year == 4)
      contour.df <- cbind(rescalee(x.df), prob = a4$prob)
  } else if (mod == "temporal") {
    current.period <- as.integer(mean.year + unq.year[year])
    current.label <- paste0("Year: ", current.period)
    alphaa <- 0.95
    contour.df <- cbind(rescalee(x.df), prob = b[[year]])
  }

  # Add first layer ------------------------------------------------------------
  p <- ggplot(contour.df, aes(x, y)) +
    geom_raster(aes(alpha = alphaa * contour.df[, 2 + 1]), fill = fill.col[1]) +
    scale_alpha_continuous(range = c(0, alphaa))

  # Add subsequent layers ------------------------------------------------------
  for (j in 2:4) {
    p <- p +
      annotate(geom = "raster", x = contour.df$x, y = contour.df$y,
               alpha = alphaa * contour.df[, 2 + j], fill = fill.col[j])
  }

  # Add decision boundaries ----------------------------------------------------
  tmp.df <- reshape2::melt(contour.df, id.vars = c("x", "y"))
  p <- p +
    geom_contour(data = tmp.df, aes(x, y, z = value, col = variable), size = 1,
                 linetype = "dashed", binwidth = 0.1,
                 breaks = seq(0.5, 0.5, by = 0.1)) +
    scale_colour_manual(values = iprior::gg_colour_hue(5)[1:4])

  # Add points -----------------------------------------------------------------
  if (isTRUE(points)) {
    if (mod == "period") {
      as.tibble(btb) %>%
        subset(sp != "Others" & period == current.period) -> points.df
    } else if (mod == "temporal") {
      as.tibble(btb) %>%
        subset(sp != "Others" & year == current.period) -> points.df
    }
    p <- p +
      geom_point(data = points.df, aes(x, y, fill = sp), col = "grey15",
                 shape = 21, size = 2) +
      scale_fill_manual(values = fill.col[1:4])
  }

  p +
    geom_polygon(data = plot.df, aes(x, y, group = id), fill = NA,
                 col = "grey25") +
    labs(x = "Eastings (1,000 km)", y = "Northings (1,000 km)") +
    scale_x_continuous(labels = function(x) x / 1000) +
    scale_y_continuous(labels = function(x) x / 1000) +
    theme_bw() +
    guides(colour = FALSE, alpha = FALSE,
           fill = guide_legend(nrow = 2,
                               title = "Spoligotype",
                               direction = "horizontal",
                               title.position = "top",
                               override.aes = list(size = 2))) +
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
          legend.title = element_text(size = 8), legend.title.align = 0.5,
          panel.grid = element_blank()) +
    # legend.box.background = element_rect(fill = NA)) +
    annotate("text", x = min(plot.df$x), y = max(plot.df$y),
             label = current.label, vjust = 1, hjust = 0)
}

p1a <- plot_stemporal_model(1)
p2a <- plot_stemporal_model(2)
p3a <- plot_stemporal_model(3)
p4a <- plot_stemporal_model(4)
cowplot::plot_grid(p1a, p2a, p3a, p4a, ncol = 2,
                   labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 10,
                   label_fontface = "plain")

# GIF
makeplot <- function() {
  for (i in seq_along(unq.year)) {
    p <- plot_stemporal_model(i, mod = "temporal")
    print(p)
  }
  animation::ani.pause()
}
animation::saveGIF(makeplot(), interval = 1, ani.width = 600, ani.height = 550)

# Bootstrap standard errors (takes a long time). Uncomment which model to fit
B <- 30
snow.options.list <- list(progress = function(i) setTxtProgressBar(pb, i))
pb <- txtProgressBar(min = 0, max = B, style = 1)
cl <- parallel::makeCluster(parallel::detectCores())
doSNOW::registerDoSNOW(cl)
res <- foreach::`%dopar%`(
  foreach::foreach(
    i = seq_len(B),
    .packages = c("iprior", "iprobit"),
    .options.snow = snow.options.list
  ), {
    ts <- sample(seq_len(length(Spoligotype)), replace = TRUE)
    # # Constant only model
    # iprobit(y = Spoligotype, X, control = list(int.only = TRUE,
    #                                            theta0 = -100))$param.full
    # Spatial model
    # iprobit(y = Spoligotype, X, kernel = "fbm",
    #         control = list(maxit = 200))$param.full
    # # Spatio-period model (period is categorical)
    # iprobit(y = Spoligotype, X, period, kernel = "fbm", interactions = "1:2",
    #         control = list(maxit = 200))$param.full
    # # Spatio-temporal model (year is continuous)
    iprobit(y = Spoligotype, X, year, kernel = "fbm", interactions = "1:2",
            control = list(maxit = 200))$param.full
  }
)
parallel::stopCluster(cl)
apply(simplify2array(res), 1:2, sd)
```