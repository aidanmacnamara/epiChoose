
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


# BLUEPRINT DATA ----------------------------------------------------------

load("r_data/blueprint_chip_cut_final.RData")
load("r_data/blueprint_rna_cut_final.RData")


# GSK DATA ----------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3","CTCF")
gsk_input = "data/data_gsk.csv"

require(BiocParallel)
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=4))

# confirm all chip matrices have the same sample order
gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)

# this will be input
 
gene_list = roi$gene
load("data/rna/E-MTAB-4101-atlasExperimentSummary.Rdata")
r_data = experimentSummary
r_data_counts = assays(r_data[[1]])$counts
r_metadata = read_tsv("data/rna/E-MTAB-4101.sdrf.txt")
chip_labels = gsk_chip_filtered[[1]]$annot$Label #
rna_labels = toupper(str_extract(r_metadata$`Source Name`, "[[:alnum:]]+_[[:alnum:]]+"))

# expression data
gsk_rna = prep_gsk_rna(r_data_counts, r_metadata, gene_list, chip_labels, rna_labels)


# INTEGRATION -------------------------------------------------------------

# TODO merge blueprint and gsk here
gsk_data = c(gsk_chip_filtered, list(gsk_rna))
blueprint_data = c(blueprint_chip_cut_final, list(blueprint_rna_cut_final))
all_data = 



pca_data = prep_for_pca_plot


# produce data frame
all_data = chip_rna_int(gsk_chip_filtered, gsk_rna, roi, marks, chip_labels, "GSK", log_data=FALSE)

# explore
dat = group_by(all_data, sample) %>% nest() # group by donor / cell type
dat_plot = dat$data[[5]] # pick first donor / cell type
ggpairs(dplyr::select_(dat_plot, .dots=names(dat_plot)[1:4])) + theme_thesis()