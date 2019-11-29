
require(tidyverse)
require(limma)
require(vsn)

# look at specific comparisons
# start with thp-1 pma vs. baseline vs. novakovic primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# NEW PRIMARY -------------------------------------------------------------

# PICK SAMPLES ------------------------------------------------------------

# export to igv to check signal noise
signal_check_ix = c(5,6,20,31,33,55,70,150,164)
dat_all$tss$H3K27ac$annot[signal_check_ix,] %>% select(`Sub-group`)

donors = list(
  monocyte = c("C0010K","C00408","C000S5","C0011I","C004SQ"), # C001UY removed
  macrophage = c("S0022I","S00390"), # S001S7 removed
  inflamm = c("S0022I","S001MJ") # S001S7 removed
)

table(unlist(donors) %in% dat_all$tss$H3K27ac$annot$Donor)

row_ix = which(
  grepl("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | (dat_all$tss$H3K27ac$annot$Donor %in% unlist(donors))
)
row_ix = c(row_ix, which(rownames(dat_all$tss$H3K27ac$res)=="Monocytes-CD14+_Broad"))

# check sample correlations
samp_corr = dat_all$tss$H3K27ac$res[(dat_all$tss$H3K27ac$annot$Donor %in% unlist(donors)), sample(1:dim(dat_all$tss$H3K27ac$res)[2], size=1e3)]
samp_corr[1:5,1:5]
rownames(samp_corr) = str_extract(rownames(samp_corr), "monocyte|macrophage|inflammatory")
samp_corr = log(samp_corr, base=2)
samp_corr[1:5,1:5]
samp_corr_df = data.frame(t(samp_corr))
ggpairs(samp_corr_df)
cor_res = cor(t(samp_corr))
corrplot(cor_res)

col_data = dat_all$tss$H3K27ac$annot[row_ix,]
col_data = data.frame(dplyr::select(col_data, Label))
col_data$rep = str_extract(col_data$Label, "BR[12]")
col_data$condition = str_extract(col_data$Label, "[[:alnum:]+]+$")
col_data$cell_type = str_extract(col_data$Label, "^[[:alnum:]-]+")
primary_ix = c(1:9,34)
col_data$rep[primary_ix] = "BR1"
col_data$cell_type[primary_ix] = "primary"
col_data$condition[c(1:5,34)] = "primary_monocyte"
col_data$condition[c(8,9)] = "primary_macrophage"
col_data$condition[c(6,7)] = "primary_macrophage_inflamm"
col_data$source = c(rep("External",9), rep("GSK",24), "External")
col_data$group = paste(col_data$cell_type, col_data$condition, sep="_")
col_data$group[primary_ix] = col_data$condition[primary_ix]
col_data[,2:6] = lapply(col_data[,2:6], factor)


# NORMALISATION -----------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix,]
dat_add[is.na(dat_add)] = 1
# dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dat_add = t(dat_add)
dat_add[1:5,1:5]

# dat_add_voomed = voom(dat_add, plot=TRUE)
# dat_add_voomed$E[1:5,1:5]

dat_add_vsn = justvsn(dat_add)
limma::plotMA(dat_add)
limma::plotMA(dat_add_vsn)
dat_add_vsn[1:5,1:5]

group = col_data$group
group = str_replace_all(group,"[-+]","_")
mm = model.matrix(~0 + group)
dds_all = lmFit(dat_add_vsn, design=mm)
rld_all = dat_add_vsn

save(dds_all, file="tmp/dds_all.RData")
save(rld_all, file="tmp/rld_all.RData")


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(rld_all))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(col_data$condition, col_data$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_all, annot_1=rld_all$condition, annot_2=rld_all$cell_type)


