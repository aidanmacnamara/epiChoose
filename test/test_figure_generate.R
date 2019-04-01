# look at specific comparisons
# start with thp-1 pma vs. baseline vs. novakovic primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# CYTOKINE DATA -----------------------------------------------------------

# cyto = read_excel("tmp/20April16 CTTV020 Project 2 MSD cytokines.xlsx")


# NOVAKOVIC ---------------------------------------------------------------

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
col_data$treatment[col_data$treatment=="Attached"] = "Naive"
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
dds_novakovic = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~treatment+donor+time)

design(dds_novakovic)
dds_novakovic$time_treatment = factor(paste(dds_novakovic$time, dds_novakovic$treatment, sep="_"))
design(dds_novakovic) = formula(~time_treatment)
dds_novakovic = DESeq(dds_novakovic)

rld_novakovic = vst(dds_novakovic, blind=FALSE)

save(dds_novakovic, file="tmp/dds_novakovic.RData")
save(rld_novakovic, file="tmp/rld_novakovic.RData")


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_novakovic)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_novakovic$treatment, rld_novakovic$time, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_novakovic, annot_1=paste(rld_novakovic$treatment, rld_novakovic$time, sep="_"), annot_2=rld_novakovic$donor)


# SAEED -------------------------------------------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(grepl("N000", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) & !is.na(dat_all$tss$H3K27ac$annot$Name))

col_data = data.frame(label=rownames(dat_all$tss$H3K27ac$res)[row_ix])
col_data$time = factor(str_extract(col_data$label, "[0-9]+(hr[s]*|days)"))
col_data$treatment = c("BG", "LPS", "Untreated", "Untreated", "BG", "LPS", "Untreated", "Untreated")
col_data$donor = str_extract(col_data$label, "^[[:alnum:]]+")
rownames(col_data) = rownames(dat_all$tss$H3K27ac$res)[row_ix]


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix,]
# dat_add = total_reg$H3K27ac$res[row_ix_filt,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_saeed = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data, design=~treatment+donor+time)

design(dds_saeed)
dds_saeed$time_treatment = factor(paste(dds_saeed$time, dds_saeed$treatment, sep="_"))
design(dds_saeed) = formula(~time_treatment)
dds_saeed = DESeq(dds_saeed)

rld_saeed = vst(dds_saeed, blind=FALSE)

save(dds_saeed, file="tmp/dds_saeed.RData")
save(rld_saeed, file="tmp/rld_saeed.RData")


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_saeed)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_saeed$treatment, rld_saeed$time, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_saeed, annot_1=paste(rld_saeed$treatment, rld_saeed$time, sep="_"), annot_2=rld_saeed$donor)


# CELL LINES --------------------------------------------------------------

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

design(dds_cell_lines)
dds_cell_lines$cell_condition = factor(paste(dds_cell_lines$cell_line,dds_cell_lines$condition,sep="_"))
design(dds_cell_lines) = formula(~cell_condition)
design(dds_cell_lines)
dds_cell_lines = DESeq(dds_cell_lines)

rld_cell_lines = vst(dds_cell_lines, blind=FALSE)

save(dds_cell_lines, file="tmp/dds_cell_lines.RData")
save(rld_cell_lines, file="tmp/rld_cell_lines.RData")


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld_cell_lines)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_cell_lines$condition, rld_cell_lines$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_cell_lines, annot_1=rld_cell_lines$condition, annot_2=rld_cell_lines$cell_line)


# ALL SAMPLES -------------------------------------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(
  grepl("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | (grepl("sanquin", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) & !is.na(dat_all$tss$H3K27ac$annot$Name)) | grepl("N000", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) 
)
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
col_data$treatment[1:8] = c("BG", "LPS", "Untreated", "Untreated", "BG", "LPS", "Untreated", "Untreated")
col_data$donor = str_replace(col_data$label, "^.*mono_([0-9]+).*$", "\\1")
col_data$donor[1:8] = str_extract(col_data$label[1:8], "^[[:alnum:]]+")
col_data$treatment[col_data$treatment=="RPMI"] = "Naive"
rownames(col_data) = rownames(dat_all$tss$H3K27ac$res)[row_ix]
col_data[length(row_ix),2:5] = c("NA",0,"Naive","Broad")

# remove challenge data
c_ix = which(col_data$challenge==1)
row_ix_filt = row_ix[-c_ix]
col_data_filt = col_data[-c_ix,]
col_data_filt$treatment[33:57] = str_extract(col_data_filt$label[33:57], "[[:alnum:]+]+$")
col_data_filt$time[33:57] = "0days"
col_data_filt$donor[33:57] = "None"
col_data_filt$cell_type = str_extract(col_data_filt$label, "macrophage|monocyte|THP-1|U937")
col_data_filt$cell_type[dim(col_data_filt)[1]] = "monocyte"


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix_filt,]
# dat_add = total_reg$H3K27ac$res[row_ix_filt,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
boxplot(t(dat_add))

dds_all = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~treatment+time)

rld_all = vst(dds_all, blind=FALSE)

save(rld_all, file="tmp/rld_all.RData")

