---
title: "California map"
date: 2018-01-28T21:55:52+01:00
anchor: "calimap"
weight: 45
---

The code to produce the map of the locations in California of interest pertaining to the Ozone data set. The latitude and longitude of the places are given in `ozone_map.txt`.

```R
# Prepare data points using ggmap package
ozone.points <- read.table("data/ozone_map.txt", header = TRUE)
ozone.points <- cbind(ozone.points, labels = c(
  "El Monte, CA",
  "Sandberg, CA",
  "Upland, CA",
  "Vandenberg AFB",
  "LAX airport",
  "Dagget, CA"
))
ozone.box <- make_bbox(lat = ozone.points$lat, lon = ozone.points$lon, f = 0.2)
ozone.box["bottom"] <- ozone.box["bottom"] - 0.2
ozone.box["top"] <- ozone.box["top"] + 0.2
ozone.map <- get_map(location = ozone.box, maptype = "watercolor",
                     crop = TRUE, source = "stamen", zoom = 10)

# Plot
ggmap(ozone.map) +
  geom_label_repel(data = ozone.points, aes(x = lon, y = lat, label = labels),
                   box.padding = 1, nudge_x = -0.05, col = "grey30",
                   segment.colour = "grey30") +
  geom_point(data = ozone.points, mapping = aes(x = lon, y = lat),
             fill = "grey30", colour = "grey30", size = 3,
             shape = 21) +
  theme_void()  # ignore warning messages
```