#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

idepth_file <- "output/020_stats/pools/stats.idepth"

idepth <- fread(idepth_file)

n_indiv <- idepth[, length(unique(INDV))]

ggplot(idepth, aes(x = MEAN_DEPTH)) +
    ylab("Number of individuals") +
    geom_histogram(bins = n_indiv)
        
