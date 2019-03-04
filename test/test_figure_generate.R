# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanquin primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# SANQUIN -----------------------------------------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(grepl("sanquin", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) & !is.na(dat_all$tss$H3K27ac$annot$Name))

col_data = data.frame(label=rownames(dat_all$tss$H3K27ac$res)[row_ix])
my_time = str_extract_all(col_data$label, "[0-9]+(hr[s]*|days)")
col_data$time = NA
col_data$challenge = 0
for(i in 1:length(my_time)) {
  if(grepl("LPS_T=4hrs$", col_data$label[i]) & length(my_time[[i]])>1) {
    col_data$challenge[i] = 1
    col_data$time[i] = my_time[[i]][2]
  } else if(length(my_time[[i]])==2) {
    col_data$time[i] = my_time[[i]][2]
  } else {
    col_data$time[i] = my_time[[i]][1]
  }
}
col_data$treatment = str_replace(col_data$label, "^.*?([[:alpha:]_]+)_T.*$", "\\1")
col_data$treatment = str_replace(col_data$treatment, "RPMI_", "")
col_data$donor = str_replace(col_data$label, "^.*mono_([0-9]+).*$", "\\1")
col_data$treatment[col_data$treatment=="RPMI"] = "Naive"
rownames(col_data) = rownames(dat_all$tss$H3K27ac$res)[row_ix]

# remove challenge data
c_ix = which(col_data$challenge==1)
row_ix_filt = row_ix[-c_ix]
col_data_filt = col_data[-c_ix,]


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix_filt,]
# dat_add = total_reg$H3K27ac$res[row_ix_filt,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_sanquin = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~treatment+donor+time)
rld_sanquin = vst(dds_sanquin, blind=FALSE)


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_sanquin)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_sanquin$treatment, rld_sanquin$time, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(filename="c:/Users/am673712/Dropbox/OTAR020/images/sanquin_image_1.png", width=1000, height=600)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)
dev.off()

png(filename="c:/Users/am673712/Dropbox/OTAR020/images/sanquin_image_2.png", width=1000, height=600)
pca_and_plot(rld_sanquin, annot_1=paste(rld_sanquin$treatment, rld_sanquin$time, sep="_"), annot_2=rld_sanquin$donor)
dev.off()




# THP-1/U937 --------------------------------------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(grepl("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE))

col_data = dat_all$tss$H3K27ac$annot[row_ix,]
col_data_filt = data.frame(dplyr::select(col_data, Label))
col_data_filt$rep = str_extract(col_data_filt$Label, "BR[12]")
col_data_filt$condition = str_extract(col_data_filt$Label, "[[:alnum:]+]+$")
col_data_filt$cell_line = str_extract(col_data_filt$Label, "^[[:alnum:]-]+")


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix,]
# dat_add = total_reg$H3K27ac$res[row_ix,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_cell_lines = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~rep+condition+cell_line)
rld_cell_lines = vst(dds_cell_lines, blind=FALSE)


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_cell_lines)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_cell_lines$condition, rld_cell_lines$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_cell_lines, annot_1=rld_cell_lines$condition, annot_2=rld_cell_lines$cell_line)



# COMPARE THP1/U937 AND SANQUIN TOGETHER ----------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(grepl("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("sanquin", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) & !is.na(dat_all$tss$H3K27ac$annot$Name))
row_ix = c(row_ix, which(rownames(dat_all$tss$H3K27ac$res)=="Monocytes-CD14+_Broad"))

col_data = data.frame(label=rownames(dat_all$tss$H3K27ac$res)[row_ix])
my_time = str_extract_all(col_data$label, "[0-9]+(hr[s]*|days)")
col_data$time = NA
col_data$challenge = 0
for(i in 1:length(my_time)) {
  if(grepl("LPS_T=4hrs$", col_data$label[i]) & length(my_time[[i]])>1) {
    col_data$challenge[i] = 1
    col_data$time[i] = my_time[[i]][2]
  } else if(length(my_time[[i]])==2) {
    col_data$time[i] = my_time[[i]][2]
  } else {
    col_data$time[i] = my_time[[i]][1]
  }
}
col_data$treatment = str_replace(col_data$label, "^.*?([[:alpha:]_]+)_T.*$", "\\1")
col_data$treatment = str_replace(col_data$treatment, "RPMI_", "")
col_data$donor = str_replace(col_data$label, "^.*mono_([0-9]+).*$", "\\1")
col_data$treatment[col_data$treatment=="RPMI"] = "Naive"
rownames(col_data) = rownames(dat_all$tss$H3K27ac$res)[row_ix]
col_data[length(row_ix),2:5] = c("NA",0,"Naive","Broad")

# remove challenge data
c_ix = which(col_data$challenge==1)
row_ix_filt = row_ix[-c_ix]
col_data_filt = col_data[-c_ix,]
col_data_filt$treatment[25:48] = str_extract(col_data_filt$label[25:48], "[[:alnum:]+]+$")
col_data_filt$time[25:48] = "0days"
col_data_filt$donor[25:48] = "None"
col_data_filt$cell_type = str_extract(col_data_filt$label, "monocyte|THP-1|U937")
col_data_filt$cell_type[dim(col_data_filt)[1]] = "monocyte"


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix_filt,]
# dat_add = total_reg$H3K27ac$res[row_ix_filt,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_all = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~treatment+time)
rld_all = vst(dds_all, blind=FALSE)


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_all)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_all$cell_type, rld_all$treatment, rld_all$time, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_all, annot_1=paste(rld_all$cell_type, rld_all$treatment, rld_all$time, sep="_"), annot_2=factor(c(rep(1,48),2)))

