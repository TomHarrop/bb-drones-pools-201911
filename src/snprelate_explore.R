library(SNPRelate)
library(data.table)
library(ggplot2)

vcf_file <- "output/020_filtered-genotypes/filtered.vcf.gz"
# vcf_file <- "test/filtered.vcf.gz"
gds_file <- tempfile(fileext = ".gds")


snpgdsVCF2GDS(vcf_file, gds_file, method = "copy.num.of.ref" )
# snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only" )
gds_data <- snpgdsOpen(gds_file)

# prune for LD and exlude non-autosomes
pruned <- snpgdsLDpruning(gds_data,
                autosome.only = FALSE, ld.threshold = 0.2)
keep_snps <- unlist(pruned[startsWith(names(pruned), "chrNC")])

# run identity by state
ibs_res <- snpgdsIBS(gds_data,
                     snp.id = keep_snps,
                     autosome.only = FALSE)

loc <- cmdscale(1 - ibs_res$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
plot(x, y, xlab = "", ylab = "", main = "cmdscale(IBS Distance)")

# convert to dist for clustering
ibs_df <- data.frame(ibs_res$ibs, row.names = ibs_res$sample.id)
colnames(ibs_df) <- ibs_res$sample.id

# cluster the ibs results
d <- as.dist(1 - ibs_df)
hc <- hclust(d, method = "ward.D2")
indiv_order <- hc$labels[hc$order]

# generic plot of the clusters
xw <- grid::convertUnit(unit(210 - 20, "mm"), "in", valueOnly = TRUE)
cairo_pdf("dendro.pdf", pointsize = 8, width = xw, height = xw)
plot(hc)
dev.off()

# data.table of ibs results for plotting
ibs_dt <- data.table(ibs_df, keep.rownames = TRUE)
setnames(ibs_dt, "rn", "i1")
ibs_pd <- melt(ibs_dt, id.vars = "i1", variable.name = "i2")

ibs_pd[, i1 := factor(i1, levels = indiv_order)]
ibs_pd[, i2 := factor(i2, levels = indiv_order)]

# set up label colours
# label_types <- gsub("[[:digit:]]+", "", indiv_order)
# label_cols <- plyr::mapvalues(label_types,
#                               from = unique(label_types),
#                               to = RColorBrewer::brewer.pal(length(unique(label_types)), "Set1"))

# plot heatmap of similarity
gp <- ggplot(ibs_pd, aes(x = i1, y = i2, fill = value)) +
    theme_grey(base_size = 8) +
    # theme(axis.text.x = element_text(angle = 90,
    #                                  hjust = 1,
    #                                  vjust = 0.5,
    #                                  colour = label_cols),
    #       axis.text.y = element_text(colour = label_cols)) +
    xlab(NULL) + ylab(NULL) + coord_fixed() +
    scale_fill_viridis_c(guide = guide_colorbar(title = "IBS"),
                         limits = c(0,1)) +
    geom_raster()

ggsave("test.pdf",gp, width = 210 - 20, height = 210 - 20, units = "mm")

# what do the GRMs look like?
grm_res <- snpgdsGRM(gds_data,
                     # method = "Weighted",
                     autosome.only = FALSE,
                     num.thread = 8,
                     snp.id = keep_snps)

# Eigen decomposition
eig <- eigen(grm_res$grm)
plot(eig$vectors[,1], eig$vectors[,2])
vect <- eig$vectors
rownames(vect) <- grm_res$sample.id
colnames(vect) <- grm_res$sample.id

# convert to dist for clustering
grm_df <- data.frame(grm_res$grm, row.names = grm_res$sample.id)
colnames(grm_df) <- grm_res$sample.id

# plot the GRM
d <- as.dist(1 - scales::rescale(as.matrix(grm_df), to = c(0,1)))
# hc2 <- hclust(as.dist(vect))
# indiv_order <- hc2$labels[hc2$order]
hc <- hclust(d, method = "ward.D2")
indiv_order <- hc$labels[hc$order]

# generic plot of the clusters
xw <- grid::convertUnit(unit(210 - 20, "mm"), "in", valueOnly = TRUE)
cairo_pdf("dendro.pdf", pointsize = 8, width = xw, height = xw)
plot(hc)
dev.off()

# data.table of ibs results for plotting
grm_dt <- data.table(grm_df, keep.rownames = TRUE)
setnames(grm_dt, "rn", "i1")
grm_pd <- melt(grm_dt, id.vars = "i1", variable.name = "i2")

grm_pd[, i1 := factor(i1, levels = indiv_order)]
grm_pd[, i2 := factor(i2, levels = indiv_order)]

# set up label colours
# label_types <- gsub("[[:digit:]]+", "", indiv_order)
# label_cols <- plyr::mapvalues(label_types,
#                               from = unique(label_types),
#                               to = RColorBrewer::brewer.pal(length(unique(label_types)), "Set1"))

# plot heatmap of similarity
gp <- ggplot(grm_pd, aes(x = i1, y = i2, fill = value)) +
    theme_grey(base_size = 8) +
    # theme(axis.text.x = element_text(angle = 90,
    #                                  hjust = 1,
    #                                  vjust = 0.5,
    #                                  colour = label_cols),
    #       axis.text.y = element_text(colour = label_cols)) +
    xlab(NULL) + ylab(NULL) + coord_fixed() +
    scale_fill_viridis_c(guide = guide_colorbar(title = "IBS")) +
    geom_raster()
gp


# SNPrelate's PCA
pca <- snpgdsPCA(gds_data, snp.id = keep_snps, autosome.only = FALSE)
pca$varprop
plot(pca)


eig <- data.table(pca$eigenvect)
setnames(eig, paste0("EV", 1:dim(eig)[[2]]))
eig[, sample := pca$sample.id]

varprop <- data.table(variable = paste0("EV", 1:length(pca$varprop)),
           varprop = pca$varprop)
varprop[, varpct := varprop * 100]           

eig_long <- melt(eig, id.vars = "sample")
eig_pd <- merge(eig_long, varprop, all.x = TRUE, all.y = FALSE)
eig_pd[, label := paste0(variable, " (", round(varpct, 1), "%)")]
eig_pd[, sample_type:= unlist(strsplit(sample, "_"))[[2]], by = sample]
setorder(eig_pd, sample_type, sample)
eig_pd[, sample := factor(sample, levels = unique(sample))]

ggplot(eig_pd[variable %in% paste0("EV", 1:10)],
       aes(x = sample, y = value, colour = sample_type)) +
    facet_wrap(~label) +
    geom_point()

