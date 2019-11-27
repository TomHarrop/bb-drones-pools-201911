#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)


frq_file <- "output/020_stats/pools/stats.frq"

frq <- fread(frq_file,
           fill = TRUE)

frq[, maf := min(`{FREQ}`,
                 V6, na.rm = TRUE),
    by = .(CHROM, POS)]

ggplot(frq, aes(x = maf)) +
    xlab("Minor allele frequency") + 
    ylab("Number of loci") + 
    geom_vline(xintercept = 0.1) +
    geom_histogram(binwidth = 0.01)
