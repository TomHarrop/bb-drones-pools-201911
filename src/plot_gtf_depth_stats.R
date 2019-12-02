#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(ggplot2)


fai_file <- "output/010_genotypes/pools/015_ref/ref.fasta.fai"
fai <- fread(fai_file)

GetNCFragments <- function(gr){
    gr[startsWith(as.character(seqnames(gr)), "NC_")]    
}

gtf_files <- list.files("output/040_phased-indivs",
                        full.names = TRUE,
                        recursive = FALSE,
                        pattern = ".gtf$")

names(gtf_files) <- sub(".gtf", "", basename(gtf_files))

gtf_list <- lapply(gtf_files, import.gff, format = "gtf")

gtf_nc <- lapply(gtf_list, GetNCFragments)

gtf_dts <- lapply(gtf_nc, function(x)
    data.table(width = width(x)))

pd <- rbindlist(gtf_dts, idcol = "indiv")
pd[, l4_width := log(width, 4)]


ggplot(pd, aes(x = l4_width)) +
    scale_x_continuous(labels = function(x) 4^x) +
    xlab("Block size") + ylab("Number of phased blocks") +
    facet_wrap(~ indiv) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = log(150, 4))

# set up bins
n_bins <- 50
weighted_pd <- pd[, .(bases = .N * width), by = .(indiv, width)]
weighted_pd[, l4_width := log(width, 4)]
quantiles <- weighted_pd[, seq(min(l4_width),
                               max(l4_width),
                               length.out = n_bins)]
bin_labels <- sapply(1:length(quantiles), function(i)
    mean(c(quantiles[i], quantiles[i + 1])))[c(1:n_bins - 1)]
weighted_pd[, wbin := as.numeric(
    as.character(
        cut(l4_width, breaks = quantiles,
            labels = bin_labels,
            include.lowest = TRUE)))]

# plot weighted histrogram
whist <- weighted_pd[, .(total_bases = sum(bases)),
                     by = .(indiv, wbin)]

ggplot(whist , aes(x = wbin, y = total_bases/1e3)) +
    scale_x_continuous(labels = function(x) 4^x) +
    facet_wrap(~ indiv) +
    ylab("Kilobases") + xlab("Haplotype block size") +
    geom_col() +
    geom_vline(xintercept = log(150, 4))

# calculate the proportion in blocks longer than x
genome_size <- fai[startsWith(V1, "NC_"), sum(V2)]

setkey(whist, wbin)
whist[, cs := cumsum(total_bases), by = .(indiv)]
whist[, genome_frac := cs / genome_size]

ggplot(whist, aes(x = wbin, y = genome_frac)) +
    ylab("Fraction of genome in haplotype block") +
    xlab("Block size") +
    scale_x_continuous(labels = function(x) 4^x) +
    facet_wrap(~ indiv) +
    geom_path() 
