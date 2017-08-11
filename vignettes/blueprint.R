
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
load_all()


# GENES -------------------------------------------------------------------

# all genes in the form of a granges object
load("r_data/column_annotation/gene_list_all.RData")


# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/column_annotation/roi.RData")


# PREPARE BLUEPRINT DATA --------------------------------------------------

# write the csv file for input
blueprint_parsed = prep_blueprint_chip(blueprint_data="data/blueprint_files.tsv", root="~/links/RD-Epigenetics-NetworkData/otar_020/BLUEPRINT/", out_file="data/blueprint_parsed.csv")

marks = c("H3K27ac","H3K4me3","H3K27me3")
blueprint_input = "data/blueprint_parsed.csv"

require(BiocParallel)
# blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
blueprint_chip = vector("list", length=length(marks))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip[[i]] = make_auc_matrix(blueprint_input, roi, marks[i], "tmp/", quantile_norm=TRUE)
}

blueprint_rna = prep_blueprint_rna(quantile_norm=TRUE)


# PLOT --------------------------------------------------------------------

plot_pca(blueprint_chip[[1]]$res, blueprint_chip[[1]]$annot$`Cell type`, "Blueprint", out_file="out_1.png")
plot_pca(blueprint_chip[[2]]$res, blueprint_chip[[2]]$annot$`Cell type`, "Blueprint", out_file="out_2.png")
plot_pca(blueprint_chip[[3]]$res, blueprint_chip[[3]]$annot$`Cell type`, "Blueprint", out_file="out_3.png")


# COMBINE CHIP/RNA --------------------------------------------------------

# filter genes
blueprint_rna$res = blueprint_rna$res[,match(roi$gene, colnames(blueprint_rna$res))]

# what donor / cell type combinations are in the rna data?
rna_labels = paste(blueprint_rna$annot$Comment.donor.ID., blueprint_rna$annot$Characteristics.cell.type., sep="_")
donor_cell = unique(rna_labels)

chip_matches = lapply(blueprint_chip, function(x) match(donor_cell, x$annot$Label))
rna_matches = match(donor_cell, rna_labels)

# filter by above indices
blueprint_chip_cut = vector("list", length(blueprint_chip))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip_cut[[i]] = list(res=blueprint_chip[[i]]$res[chip_matches[[i]],], annot=blueprint_chip[[i]]$annot[chip_matches[[i]],])  
}

blueprint_rna_cut = lapply(blueprint_rna, function(x) x[rna_matches,])

# filter missing data
final_ix = unique(unlist(lapply(blueprint_chip_cut, function(x) which(apply(x[[1]], 1, function(x) !all(is.na(x)))))))
blueprint_chip_cut_final = vector("list", length(blueprint_chip))
for(i in 1:length(blueprint_chip_cut)) {
  blueprint_chip_cut_final[[i]] = list(res=blueprint_chip_cut[[i]]$res[final_ix,], annot=blueprint_chip_cut[[i]]$annot[final_ix,])  
}

blueprint_rna_cut_final = lapply(blueprint_rna, function(x) x[final_ix,])

dim_1 = dim(blueprint_rna_cut_final$res)[1]
dim_2 = dim(blueprint_rna_cut_final$res)[2]

all_data = data.frame(
  rna_response = log(as.numeric(t(blueprint_rna_cut_final$res))),
  chip_k27ac = log(as.numeric(t(blueprint_chip_cut_final[[1]]$res))),
  chip_k4me3 = log(as.numeric(t(blueprint_chip_cut_final[[2]]$res))),
  chip_k27me3 = log(as.numeric(t(blueprint_chip_cut_final[[3]]$res))),
  group = factor(rep(blueprint_chip_cut_final[[1]]$annot$Group, each=dim_2)),
  donor = factor(rep(blueprint_chip_cut_final[[1]]$annot$Donor, each=dim_2)),
  cell_type = factor(rep(blueprint_chip_cut_final[[1]]$annot$`Cell type`, each=dim_2)),
  gene = factor(rep(roi$gene, dim_1))
)

all_data[all_data==-Inf] = 0 # reconvert to 0
dat = group_by(all_data, donor, cell_type, group) %>% nest() # group by donor / cell type
dat_2 = dat$data[[33]] # pick first donor / cell type
png(filename="out.png", height=800, width=1200)
ggpairs(dplyr::select(dat_2, rna_response, chip_k27ac, chip_k4me3, chip_k27me3)) + theme_thesis()
dev.off()


