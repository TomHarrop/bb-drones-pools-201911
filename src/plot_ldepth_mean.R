#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

ldepth_mean_file <- "output/020_stats/drones/stats.ldepth.mean"

ldepth_mean <- fread(ldepth_mean_file)

mean_depth <- ldepth_mean[, mean(MEAN_DEPTH)]

# fwrite this and use for filtering
data.table(c("min_depth", "max_depth"),
           c(0.5, 2) * mean_depth)


ggplot(ldepth_mean, aes(x = MEAN_DEPTH)) +
    ylab("Number of loci") +
    scale_x_log10() +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = c(0.5, 2) * mean_depth)
