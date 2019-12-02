#!/usr/bin/env Rscript

library(VariantAnnotation)
library(data.table)
library(ggplot2)

RetrieveGtDataTable <- function(vcf_file, tbi_file) {
  my_tabix <- TabixFile(vcf_file, tbi_file)
  my_vcf <- readVcf(my_tabix)
  my_vcf_dt <- as.data.table(geno(my_vcf)$GT, keep.rownames = TRUE)
  setnames(my_vcf_dt, "rn", "locus")
  return(melt(my_vcf_dt,
              id.vars = "locus",
              variable.name = "indiv",
              value.name = "GT"))}

dr_vcf <- RetrieveGtDataTable(
  "output/020_filtered-genotypes/drones/filtered.vcf.gz",
  "output/020_filtered-genotypes/drones/filtered.vcf.gz.tbi")

po_vcf <- RetrieveGtDataTable(
  "output/020_filtered-genotypes/pools/filtered.vcf.gz",
  "output/020_filtered-genotypes/pools/filtered.vcf.gz.tbi")

merged_gt <- merge(dr_vcf, po_vcf,
                   all = TRUE,
                   by = c("locus", "indiv"),
                   suffixes = c("_dr", "_po"))

# get an overview (not by indiv)
merged_gt[, table(GT_po, GT_dr)]

# decide if loci are [homozygous] variants
merged_gt[, `drone-var` := as.numeric(GT_dr) > 0 ]
merged_gt[, `pool-hom-var` :=
            unlist(lapply(strsplit(GT_po, "/"), function(x)
              all(as.numeric(x) > 0)))]

merged_gt[, table(`drone-var`, `pool-hom-var`)]

# tests
merged_gt[GT_dr == "0" & `pool-hom-var` == TRUE,
          gt_match := "drone-ref_pool-hom-var"]
merged_gt[`drone-var` == TRUE & GT_po == "0/0",
          gt_match := "drone-var_pool-hom-ref"]
merged_gt[is.na(GT_dr) & `pool-hom-var` == TRUE,
          gt_match := "drone-na_pool-hom-var"]
merged_gt[is.na(GT_dr) & is.na(GT_po)]          # no such thing
merged_gt[`drone-var` == TRUE & is.na(`pool-hom-var`),
          gt_match := "drone-var_pool-na"]
merged_gt[is.na(gt_match),
          gt_match := "compatible"]

merged_gt[is.na(gt_match), table(GT_po, GT_dr, useNA = "always")]
merged_gt[gt_match == "drone-var_pool-na"]

merged_gt[gt_match == "match", table(GT_dr, GT_po, useNA = "always")]




# summarise
match_dt <- merged_gt[, .(n_loci = length(unique(locus))),
          by = .(indiv, gt_match)]
ggplot(match_dt, aes(x = gt_match, y = n_loci, colour = indiv)) +
  scale_color_viridis_d() +
  scale_y_log10() +
  geom_point(position = "jitter")








