
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
load_all()


# ROI ---------------------------------------------------------------------

# look across all regulatory regions
load("r_data/roi_ensembl_multicell.RData")

marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC")


# BLUEPRINT DATA ----------------------------------------------------------

blueprint_input = "data/data_blueprint_parsed.csv"

require(BiocParallel)
# blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
blueprint_chip = vector("list", length=length(marks))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip[[i]] = make_auc_matrix(blueprint_input, roi, marks[i], "tmp/", quantile_norm=TRUE)
}

# make sure rows are in the same order
blueprint_chip_filtered = prep_gsk_chip_filter(blueprint_chip)

# take out samples with missing data for 1 or more data type
na_df = data.frame(lapply(blueprint_chip_filtered, function(x) apply(x$res, 1, function(y) !all(is.na(y)))))
names(na_df) = NULL
na_ix = which(apply(na_df, 1, all))

for(i in 1:length(blueprint_chip_filtered)) {
  blueprint_chip_filtered[[i]]$res = blueprint_chip_filtered[[i]]$res[na_ix,]
  blueprint_chip_filtered[[i]]$annot = blueprint_chip_filtered[[i]]$annot[na_ix,]
}


# GSK DATA ----------------------------------------------------------------

gsk_input = "data/data_gsk.csv"
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=4))
gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)


# ENCODE DATA -------------------------------------------------------------

encode_input = "data/data_encode.csv"
encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
encode_chip_filtered = prep_gsk_chip_filter(encode_chip)


# COMBINE -----------------------------------------------------------------

mask_data = vector("list", 3)
for(i in 1:length(mask_data)) {
  mask_data[[i]]$res = rbind(blueprint_chip_filtered[[i]]$res, gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res)
  # renormalize as we are multiple sources
  mask_data[[i]]$res = quantile_norm(mask_data[[i]]$res)
  mask_data[[i]]$annot = bind_rows(blueprint_chip_filtered[[i]]$annot, gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot)
}

names(mask_data) = marks


# ANALYSIS ----------------------------------------------------------------

group_labels = c(rep("BLUEPRINT",76), rep("GSK",18), rep("ENCODE",18))
single_labels = rownames(mask_data[[2]]$res)

# total plot
pca_data = prep_for_plot(mask_data, annot_1=group_labels, annot_2=single_labels, marks=names(mask_data), plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + facet_wrap(~mark, nrow=1, scales="free") # + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5)
dev.off()

# masking based on all data types
sig_list = vector("list", length(mask_data))
names(sig_list) = names(mask_data)

for(i in 1:length(sig_list)) {
  limma_res = limma_diff(mask_data[[i]]$res, groups=list(1:16, 17:50), thresh=5000)
  limma_res$gene = roi$gene[as.numeric(rownames(limma_res))]
  limma_gsea = arrange(limma_res, desc(diff)) %>% filter(adj.P.Val<0.0001, abs(diff)>1) %>% select(gene, diff)
  write_tsv(limma_gsea, paste0("c:/Downloads/tmp/gsea/gene_ranks/otar_", i, ".rnk"), col_names=FALSE)  
  sig_list[[i]] = limma_gsea
}

venn.diagram(lapply(sig_list, function(x) x$gene), filename="out.png", imagetype="png")
venn_res = calculate.overlap(lapply(sig_list, function(x) x$gene))
venn_genes = c(venn_res$a6, venn_res$a7, venn_res$a5, venn_res$a2)
feature_ix = match(venn_genes, roi$gene)

mask_data_sliced = mask_data
for(i in 1:length(mask_data_sliced)) {
  mask_data_sliced[[i]]$res = mask_data_sliced[[i]]$res[,-feature_ix]
}

pca_data = prep_for_plot(mask_data_sliced, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

