
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
load_all()


# GENES -------------------------------------------------------------------

# all genes in the form of a granges object
load("r_data/gene_list_all.RData")


# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/roi.RData")


# BLUEPRINT DATA ----------------------------------------------------------

load("r_data/blueprint_chip_cut_final.RData")
load("r_data/blueprint_rna_cut_final.RData")


# GSK DATA ----------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3")
gsk_input = "data/data_gsk.csv"

require(BiocParallel)
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=3))

# TO ADD confirm all chip matrices have the same sample order

prep_gsk_rna()

