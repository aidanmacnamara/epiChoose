
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
load_all()


# ALL DATA ----------------------------------------------------------------

load("r_data/all_data.RData")
load("r_data/gsk_chip.RData")
names(all_data) = c("RNA","H3k27ac","H3K4me3","H3K27me3")

# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/roi.RData")
load("r_data/gene_list_go.RData")


# ANALYSIS ----------------------------------------------------------------

gsk_chip_filtered = prep_gsk_chip_filter(gsk_chip)


# SLIDES ------------------------------------------------------------------

# 1. compare replicates

group_labels = str_replace(rownames(gsk_chip_filtered[[1]]$res), "_BR[12]", "")
single_labels = rownames(gsk_chip_filtered[[1]]$res)

plot_data = prep_for_plot(gsk_chip_filtered[1:3], annot_1=group_labels, annot_2=single_labels, marks=c("H3K27ac","H3K4me3","H3K27me3"), plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(plot_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=8, shape=17) + theme_thesis(60) + geom_text_repel(aes(label=annot_2), fontface="bold", size=10, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# 2. distance table of cell lines vs. primary cells

group_labels = c(rep("GSK",16), rep("BLUEPRINT",34))
single_labels = rownames(all_data[[2]]$res)

sample_ix = c(5:16,23,42)
all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,]
  all_data_sliced[[i]]$annot = all_data_sliced[[i]]$annot[sample_ix,]
}

plot_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"), plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(plot_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis(40) + geom_text_repel(aes(label=annot_2), fontface="bold", size=10, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

dists = vector("list", length(all_data_sliced))
names(dists) = names(all_data_sliced)

for(i in 1:length(all_data_sliced)) {
  dist_dat = filter(plot_data, mark==names(all_data_sliced)[i]) %>% dplyr::select(x, y, annot_2)
  rownames(dist_dat) = dist_dat$annot_2
  dist_dat = dist_dat[,-3]
  dists[[i]] = dist(dist_dat)
  write.csv(as.matrix(dists[[i]]), file=paste0("out_", i, ".csv"))
}

# 3. use haematopietic gene set

feature_ix = match(gene_list_go$hgnc_symbol, roi$gene)

all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[sample_ix,feature_ix]
}

lapply(all_data_sliced, function(x) dim(x[[1]]))

plot_data = prep_for_plot(all_data_sliced, annot_1=group_labels[sample_ix], annot_2=single_labels[sample_ix], marks=names(all_data_sliced), plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(plot_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

dists = vector("list", length(all_data_sliced))
names(dists) = names(all_data_sliced)

for(i in 1:length(all_data_sliced)) {
  dist_dat = filter(plot_data, mark==names(all_data_sliced)[i]) %>% dplyr::select(x, y, annot_2)
  rownames(dist_dat) = dist_dat$annot_2
  dist_dat = dist_dat[,-3]
  dists[[i]] = dist(dist_dat)
  write.csv(as.matrix(dists[[i]]), file=paste0("out_", i, ".csv"))
}
