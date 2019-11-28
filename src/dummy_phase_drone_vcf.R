library(VariantAnnotation)
library(data.table)

vcf_file <- "output/020_filtered-genotypes/drones/filtered.vcf.gz"
tbi_file <- "output/020_filtered-genotypes/drones/filtered.vcf.gz.tbi"

# read the un-phased VCF
tabix_file <- TabixFile(vcf_file, tbi_file)
vcf <- readVcf(tabix_file)

# get the genotype matrix
gt <- geno(vcf)$GT

# use RLE to make each chr a phase set (not sure if necessary)
# unlist(gt)

# generate a PS matrix using dims from genotype matrix (maybe fill with rle-encoded phase numbers)
ps <- matrix(1,
       nrow = nrow(gt),
       ncol = ncol(gt))
colnames(ps) <- colnames(gt)
rownames(ps) <- rownames(gt)

# add the PS row to the metadata
geno <- data.table(data.frame(geno(header(vcf))),
                     keep.rownames = TRUE)
ps_row <- data.table("rn" = "PS",
  "Number" = "1",
  "Type" = "Integer",
  "Description" = "Phase set")
new_geno <- data.frame(rbind(geno, ps_row), row.names = "rn")

# add the new metadata back to the vcf header
new_header <- VCFHeader(reference = reference(header(vcf)),
          samples = samples(header(vcf)))
geno(new_header) <- DataFrame(new_geno)
info(new_header) <- info(header(vcf))
fixed(new_header) <- fixed(header(vcf))

header(vcf) <- new_header

# add the PS field
geno(vcf)$PS <- ps

# write out
writeVcf(vcf,
         "test/phased_drones.vcf")

