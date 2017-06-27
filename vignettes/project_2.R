
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
load("r_data/roi.RData")


# GSK DATA ----------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3","CTCF")
gsk_input = "data/data_gsk.csv"

require(BiocParallel)
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=4))

# confirm all chip matrices have the same sample order
gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)

single_labels = gsk_chip_filtered[[2]]$annot$Label
group_labels = str_extract(single_labels, "^[[:alnum:]]+")

# total plot
pca_data = prep_for_plot(gsk_chip_filtered, annot_1=group_labels, annot_2=single_labels, marks=marks)

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()



