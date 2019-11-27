#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

qual_file <- "output/020_stats/stats.lqual"

qual <- fread(qual_file)
ggplot(qual, aes(x = QUAL)) +
  ylab("Number of loci") +
  scale_x_log10() +
  geom_vline(xintercept = 30) +
  geom_histogram(bins = 50)
