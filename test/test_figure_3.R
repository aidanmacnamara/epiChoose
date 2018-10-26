# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanquin primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# LOOK AT BLUEPRINT DATA ALONE: H3K27AC -----------------------------------

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


# SLICE AND RUN DDS/RLOG --------------------------------------------------

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

y = t(assays(rld_sanquin)[[1]])
dim(y)
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(rld_sanquin$treatment, rld_sanquin$time, sep="_")
yy$annot_2 = rld_sanquin$donor

png(filename="c:/Users/am673712/Dropbox/OTAR020/images/sanquin_image_2.png", width=1000, height=600)
ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")
dev.off()


# DIFFERENTIAL ANALYSIS ---------------------------------------------------

# novakovic criteria
# fc 3
# adj p 0.05
# rpkm > 1

design(dds_sanquin)
dds_sanquin = DESeq(dds_sanquin, test="LRT", reduced=~donor)
res_sanquin = tbl_df(results(dds_sanquin))
res_sanquin$gene = rownames(res_sanquin)
res_sanquin$symbol = gene_list_all$hgnc_symbol
arrange(res_sanquin, padj) %>% filter(padj<1e-8)

# take the top 1000 hits
genes_top = head(res_sanquin$gene[order(res_sanquin$padj)],1e3)

for(j in 1:3) {
  y = plotCounts(dds_sanquin, gene=head(order(res_sanquin$padj),10)[j], intgroup=c("treatment","time","donor"), returnData=TRUE)
  print(ggplot(y, aes(x=time, y=count)) + geom_point(shape=17, size=3) + theme_thesis(20) + xlab("") + ylab("Count") + ggtitle(res_sanquin$gene[head(order(res_sanquin$padj),10)[j]]))
}

y = t(assays(rld_sanquin)[[1]])
dim(y)
y = y[,colnames(y) %in% genes_top]
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(rld_sanquin$treatment, rld_sanquin$time, sep="_")
yy$annot_2 = rld_sanquin$donor
ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")

for_heatmap = assays(rld_sanquin)[[1]]
colnames(for_heatmap) = paste(col_data_filt$treatment, col_data_filt$time, col_data_filt$donor, sep="_")
for_heatmap = as.data.frame(for_heatmap[match(genes_top,rownames(dds_sanquin)),])
dim(for_heatmap)
pheatmap(for_heatmap, show_rownames=FALSE)

wss = (nrow(for_heatmap)-1) * sum(apply(for_heatmap,2,var))
for(k in 2:15) {
  wss[k] <- sum(kmeans(for_heatmap, centers=k)$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within-group sum-of-squares")
fit <- kmeans(for_heatmap, 6)

pheatmap(for_heatmap, fontsize_row=8,
         annotation_row = data.frame(
           K_Means = factor(fit$cluster),
           row.names=rownames(for_heatmap)
         ), show_rownames=FALSE
)

for_heatmap$Gene = rownames(for_heatmap)
for_heatmap$Cluster = factor(fit$cluster)
head(for_heatmap)
for_heatmap = gather(for_heatmap, "Group", "AUC", 1:17)
for_heatmap$Group = factor(for_heatmap$Group)
ggplot(for_heatmap, aes(x=Group,y=AUC,group=Gene)) + geom_point(shape=17) + geom_line(alpha=0.1) + facet_wrap(~Cluster) + theme_thesis(15)


# MONOCYTE/MACROPHAGE RELEVANT GENES --------------------------------------

# pick the gene set from saeed/novakovic
p_genes = read_excel("tmp/1-s2.0-S0092867416313162-mmc2.xlsx", sheet="Table S2B. Genes H3K27ac prom", skip=1)
dg_genes = p_genes$Promoter_H3K27ac_Differentiation_gain; dg_genes = dg_genes[!is.na(dg_genes)]
dl_genes = p_genes$Promoter_H3K27ac_Differentiation_loss; dl_genes = dl_genes[!is.na(dl_genes)]

for_heatmap = assays(rld_sanquin)[[1]]
colnames(for_heatmap) = paste(col_data_filt$treatment, col_data_filt$time, col_data_filt$donor, sep="_")
match_ix = match(c(dg_genes,dl_genes), gene_list_all$hgnc_symbol); match_ix = match_ix[!is.na(match_ix)]
for_heatmap = as.data.frame(for_heatmap[match_ix,])
dim(for_heatmap)

pheatmap(for_heatmap, cluster_rows=FALSE, show_rownames=FALSE)
col_clust = hclust(dist(t(for_heatmap)))
col_names = paste0("Group_", cutree(col_clust, k=2))
names(for_heatmap) = paste(names(for_heatmap), col_names, sep="_")
head(for_heatmap)
for_heatmap$Gene = gene_list_all$hgnc_symbol[match_ix]
for_heatmap$Cluster = factor(ifelse(for_heatmap$Gene %in% dg_genes, "Up", "Down"))
head(for_heatmap)

for_heatmap_long = gather(for_heatmap,"Group","AUC",1:24)
for_heatmap_long$Group = factor(for_heatmap_long$Group)
head(for_heatmap_long)
for_heatmap_long$all_group = str_extract(for_heatmap_long$Group, "Group_[12]")
ggplot(for_heatmap_long, aes(x=all_group,y=AUC)) + geom_boxplot() + facet_wrap(~Cluster) + theme_thesis(20)


# BRING IN THP-1 PMA VS. BASELINE -----------------------------------------

# pick samples
row_ix = which(grepl("thp-1", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE))

col_data = dat_all$tss$H3K27ac$annot[row_ix,]
col_data_filt = data.frame(dplyr::select(col_data, Label))
col_data_filt$rep = str_extract(col_data_filt$Label, "BR[12]")
col_data_filt$condition = str_extract(col_data_filt$Label, "[[:alnum:]+]+$")

# run dds
dat_add = dat_all$tss$H3K27ac$res[row_ix,]
# dat_add = total_reg$H3K27ac$res[row_ix,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_thp = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~rep+condition)
rld_thp = vst(dds_thp, blind=FALSE)

# plot
sampleDists <- dist(t(assay(rld_thp)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_thp$condition, rld_thp$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

y = t(assays(rld_thp)[[1]])
dim(y)
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = rld_thp$condition

ggplot(yy, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis()

design(dds_thp)
dds_thp = DESeq(dds_thp)
res_thp = tbl_df(results(dds_thp, contrast=c("condition","PMA","Baseline")))
res_thp$gene = rownames(res_thp)
res_thp$symbol = gene_list_all$hgnc_symbol
arrange(res_thp, padj)


# SHOW THE OVERALL DISTANCE -----------------------------------------------

# pick samples
row_ix = which(grepl("thp-1", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("sanquin", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) & !is.na(dat_all$tss$H3K27ac$annot$Name))

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
col_data_filt$treatment[25:36] = str_extract(col_data_filt$label[25:36], "[[:alnum:]+]+$")
col_data_filt$time[25:36] = "0days"
col_data_filt$donor[25:36] = "None"

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
rownames(sampleDistMatrix) = paste(rld_all$treatment, rld_all$time, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

y = t(assays(rld_all)[[1]])
dim(y)
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(rld_all$treatment, rld_all$time, sep="_")
yy$annot_2 = rld_all$donor

ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")

# look across go terms for highest correlation
design(dds_sanquin)
dds_sanquin$time_treatment = factor(paste(dds_sanquin$time, dds_sanquin$treatment, sep="_"))
design(dds_sanquin) = formula(~time_treatment+donor)
dds_sanquin = DESeq(dds_sanquin)
res_primary_diff = results(dds_sanquin, contrast=c("time_treatment","6days_Naive","1hr_Attached"))
res_celline_diff = results(dds_thp, contrast=c("condition","PMA","Baseline"))
dat = data.frame(prim_fc=res_primary_diff$log2FoldChange, cell_fc=res_celline_diff$log2FoldChange, gene=gene_list_all$hgnc_symbol)

# get go lists
load("../epiView/data/msig_go_bp.RData")
go_res = data.frame(cor=rep(NA, length(msig_go_bp)), name=names(msig_go_bp))

for(i in 1:dim(go_res)[1]) {
  g_ix = which(dat$gene %in% msig_go_bp[[i]])
  if(length(g_ix)<20) next
  go_res$cor[i] = cor.test(dat$prim_fc[g_ix], dat$cell_fc[g_ix])$estimate  
}

arrange(go_res, desc(cor)) %>% DT::datatable()

g_ix = which(gene_list_all$hgnc_symbol %in% msig_go_bp[[which(names(msig_go_bp)=="REGULATION OF TYPE 2 IMMUNE RESPONSE")]])

y = t(assays(dds_all)[[1]])
y = y[,g_ix]
dim(y)
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(rld_all$treatment, rld_all$time, sep="_")
yy$annot_2 = rld_all$donor

ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")

y = assays(dds_all)[[1]]
y = data.frame(y[g_ix,])
names(y) = paste(rld_all$treatment, rld_all$time, sep="_")
y$gene = gene_list_all$hgnc_symbol[g_ix]
y = gather(data.frame(y),"Group","AUC", 1:36)
y = filter(y, Group %in% c("Naive_6days","Naive_6days.1","Attached_1hr","Attached_1hr.1"))

ggplot(y, aes(x=Group,y=AUC)) + geom_point(shape=17) + geom_line(alpha=0.1) + facet_wrap(~gene) + theme_thesis(15)
