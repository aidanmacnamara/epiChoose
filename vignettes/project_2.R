
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
gsk_data = c(list(gsk_rna), gsk_chip_filtered)
blueprint_data = c(list(blueprint_rna_cut_final), blueprint_chip_cut_final)
all_data = vector("list", 4)
for(i in 1:4) {
  all_data[[i]]$res = rbind(gsk_data[[i]]$res, blueprint_data[[i]]$res)
  # renormalize as we are combining gsk and blueprint data
  all_data[[i]]$res = quantile_norm(all_data[[i]]$res)
  all_data[[i]]$annot = bind_rows(gsk_data[[i]]$annot, blueprint_data[[i]]$annot)
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


# DATA SLICES -------------------------------------------------------------

# start off with all_data, but slice by gene list / samples of interest

# 2. slice by relevant samples and whole genome
sample_ix = c(5:16,23,42)
all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,]
  all_data_sliced[[i]]$annot = all_data_sliced[[i]]$annot[sample_ix,]
}

pca_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# 3/4. mask cell line specific loadings
feature_ix = limma_diff(all_data[[2]]$res, groups=list(1:16, 17:50), thresh=5000)
# feature_ix = pca_drivers(all_data[[2]]$res, thresh=5000, which_pc=1)
all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[,-feature_ix]
}

pca_data = prep_for_plot(all_data_sliced, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# 5. slice by relevant samples and expanded gene set
# use gene list from 'CTTV/examples/poc/poc_1.R'
load("r_data/gene_list_go.RData")
feature_ix = match(gene_list_go$hgnc_symbol, roi$gene)

all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,feature_ix]
}

pca_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# 6. slice by relevant samples and "erythroid differentiation"

genes = gene_list_go$hgnc_symbol[gene_list_go$go_id %in% c("GO:0045646","GO:0030218","GO:0060319","GO:0043249","GO:0048821")]
feature_ix = match(genes, roi$gene)

all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,feature_ix]
}

pca_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# 7. slice by relevant samples and go terms

genes = gene_list_go$hgnc_symbol[gene_list_go$go_id %in% c("GO:0045646","GO:0030218","GO:0060319","GO:0043249","GO:0048821")]

feature_ix = match(genes, roi$gene)

all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,feature_ix]
}

pca_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()


# REDO THIS ---------------------------------------------------------------

all_data = chip_rna_int(gsk_chip_filtered, gsk_rna, roi, marks, chip_labels, "GSK", log_data=FALSE)

# explore
dat = group_by(all_data, sample) %>% nest() # group by donor / cell type
dat_plot = dat$data[[5]] # pick first donor / cell type
ggpairs(dplyr::select_(dat_plot, .dots=names(dat_plot)[1:4])) + theme_thesis()

