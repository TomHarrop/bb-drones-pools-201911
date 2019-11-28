#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

library(data.table)

vcf_file <- args[[1]]

vcf <- fread(cmd = paste('grep -v "^#"', vcf_file))

# define the columns
md_cols <- paste0("V", 1:9)
gt_cols <- names(vcf)[!names(vcf) %in% md_cols]

# add the format info
vcf[, V9 := paste(V9, "PS", sep = ":")]

# add the phase set
myfunc <- function(x) {
    paste(x, "100", sep = ":")
}
new_gt <- vcf[, lapply(.SD, myfunc), .SDcols = gt_cols]

# combine the new genotype columns with the md columns
new_vcf <- cbind(vcf[, md_cols, with = FALSE], new_gt)

# write the output
fwrite(new_vcf,
       args[[2]],
       sep = "\t",
       col.names = FALSE)

# log
sessionInfo()
