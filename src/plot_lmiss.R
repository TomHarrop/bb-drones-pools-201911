#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

lmiss_file <- "output/020_stats/pools/stats.lmiss"

lmiss <- fread(lmiss_file)

ggplot(lmiss, aes(x = F_MISS)) +
    xlab("Number of loci") +
    geom_histogram(bins=50) +
    geom_vline(xintercept = 0.1)
