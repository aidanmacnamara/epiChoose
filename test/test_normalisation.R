
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(BiocParallel)
require(SummarizedExperiment)
load_all()


# GENES -------------------------------------------------------------------

# all genes in the form of a granges object
load("r_data/gene_list_all.RData")


# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/roi.RData")


# COMPARE UNNORMALISED ----------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3")

gsk_input = "data/data_gsk.csv"
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=4))
gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)

blueprint_parsed = prep_blueprint_chip(blueprint_data="data/blueprint_files.tsv", root="~/links/RD-Epigenetics-NetworkData/otar_020/BLUEPRINT/", out_file="data/data_blueprint_parsed.csv")
blueprint_input = "data/data_blueprint_parsed.csv"
blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))

encode_input = "data/data_encode.csv"
encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))


# INTEGRATION -------------------------------------------------------------

# merge data
all_data = vector("list", 4)
for(i in 1:3) {
  all_data[[i]]$res = rbind(gsk_chip[[i]]$res, blueprint_chip[[i]]$res, encode_chip[[i]]$res)
  # renormalize
  # all_data[[i]]$res = quantile_norm(all_data[[i]]$res)
  all_data[[i]]$annot = bind_rows(gsk_chip[[i]]$annot, blueprint_chip[[i]]$annot, encode_chip[[i]]$annot)
}

group_labels = c(rep("GSK",16), rep("BLUEPRINT",34))
single_labels = rownames(all_data[[2]]$res)

# test plot
plot_pca(all_data[[2]]$res, annot_1=group_labels, annot_2=rownames(all_data[[2]]$res), out_file="out.png")

# total plot
pca_data = prep_for_plot(all_data, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

