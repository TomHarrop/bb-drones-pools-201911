#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

imiss_file <- "output/020_stats/pools/stats.imiss"

imiss <- fread(imiss_file)

ggplot(imiss, aes(x = F_MISS)) +
    geom_histogram(binwidth = 0.01)
