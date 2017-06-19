
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
load_all()


# GENES -------------------------------------------------------------------

# all genes in the form of a granges object
load("r_data/gene_list_all.RData")


# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/roi_ensembl_multicell.RData")


# GSK DATA ----------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3","CTCF")
gsk_input = "data/project_3.csv"

require(BiocParallel)
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=4))

# confirm all chip matrices have the same sample order
gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)

# add encode a549
load("r_data/test_masking/encode_chip_filtered.RData")

all_data = vector("list", 3)
for(i in 1:length(all_data)) {
  all_data[[i]]$res = rbind(gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res[1:3,])
  # renormalize as we are multiple sources
  all_data[[i]]$res = quantile_norm(all_data[[i]]$res)
  all_data[[i]]$annot = bind_rows(gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot[1:3,])
}

group_labels = c(str_replace(rownames(gsk_chip_filtered[[1]]$res), "_BR[12]", ""), rep("ENCODE",3))
single_labels = rownames(all_data[[1]]$res)

plot_data = prep_for_plot(all_data, annot_1=group_labels, annot_2=single_labels, marks=marks, plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(plot_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=8, shape=17) + theme_thesis(60) + geom_text_repel(aes(label=annot_2), fontface="bold", size=10, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

