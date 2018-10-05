# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanquin primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# LOOK AT BLUEPRINT DATA ALONE: H3K27AC -----------------------------------

# PICK SAMPLES ------------------------------------------------------------

row_ix = which((grepl("monocyte", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("macrophage", rownames(dat_all$tss$H3K27ac$res))) & !is.na(dat_all$tss$H3K27ac$annot$Label))

col_data = dat_all$tss$H3K27ac$annot[row_ix,]

col_data_filt = data.frame(select(col_data, Label))
col_data_filt$cell = tolower(str_extract(col_data_filt$Label, "[Mm]onocyte|macrophage"))
col_data_filt$donor = str_extract(col_data_filt$Label, "^[CS][[:alnum:]]+")
col_data_filt$donor[grepl("SANQUIN", col_data$Name)] = str_extract(col_data_filt$Label[grepl("SANQUIN", col_data$Name)], "mono_[0-9]+")
col_data_filt$condition = tolower(str_extract(col_data_filt$Label, "inflammatory|classical|alternatively activated|BR.*$"))
col_data_filt$condition = str_replace(col_data_filt$condition, "br[12]_", "")
col_data_filt$condition[is.na(col_data_filt$condition)] = "primary"
col_data_filt$condition[grepl("SANQUIN", col_data$Name)] = tolower(str_extract(col_data_filt$Label[grepl("SANQUIN", col_data$Name)], "[[:alnum:]=_]+$"))

rownames(col_data_filt) = rownames(dat_all$tss$H3K27ac$res)[row_ix]


# SLICE AND RUN DDS/RLOG --------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix,]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~condition)
rld = vst(dds, blind=FALSE)


# PLOT HEATMAP/PCA --------------------------------------------------------

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = rld$Label
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

y = t(assays(rld)[[1]])
dim(y)
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(colData(rld)$cell, colData(rld)$condition, sep="_")
yy$annot_2 = colData(rld)$condition
ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=2, force=0.5) + theme(legend.position="none")


# DIFFERENTIAL ANALYSIS ---------------------------------------------------

design(dds)
dds = DESeq(dds, test="LRT", reduced=~1)
res = tbl_df(results(dds))
res$gene = rownames(res)
res$symbol = gene_list_all$hgnc_symbol
res

# take the top 1000 hits
genes_top = head(res$gene[order(res$padj)],1e3)

for_heatmap = assays(rld)[[1]]
colnames(for_heatmap) = paste(col_data_filt$cell, col_data_filt$condition, sep="_")
for_heatmap = as.data.frame(for_heatmap[match(genes_top,rownames(dds)),])
dim(for_heatmap)
pheatmap(for_heatmap)
for_heatmap <- data.frame(
  sapply(unique(names(for_heatmap)), # for each unique column name
         function(col) rowMeans(for_heatmap[names(for_heatmap)==col]) # calculate row means
  )
)
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

for_heatmap = assays(rld)[[1]]
colnames(for_heatmap) = col_data_filt$Label
match_ix = match(c(dg_genes,dl_genes), gene_list_all$hgnc_symbol); match_ix = match_ix[!is.na(match_ix)]

for_heatmap = as.data.frame(for_heatmap[match_ix,])
dim(for_heatmap)
pheatmap(for_heatmap, cluster_rows=FALSE, show_rownames=FALSE)
fit <- kmeans(for_heatmap, 2)
pheatmap(for_heatmap, cluster_rows=FALSE, show_rownames=FALSE,
         annotation_row = data.frame(K_Means=factor(fit$cluster), row.names=rownames(for_heatmap))
)

for_heatmap$Gene = rownames(for_heatmap)
for_heatmap$Cluster = factor(ifelse(for_heatmap$Gene %in% dg_genes, "Up", "Down"))
head(for_heatmap)

# names(for_heatmap)[1:16] = paste(names(for_heatmap)[1:16],"_", 1:16, sep="")
# for_heatmap = gather(for_heatmap, "Group","AUC",1:16)
# for_heatmap$Group = factor(for_heatmap$Group)
# head(for_heatmap)
# for_heatmap$Group = str_replace(for_heatmap$Group, "_[0-9]+$", "")
# ggplot(for_heatmap, aes(x=Group,y=AUC)) + geom_boxplot() + facet_wrap(~Cluster) + theme_thesis(20)

for_heatmap$macrophage_mean = apply(for_heatmap[,names(for_heatmap)=="macrophage"], 1, mean)
for_heatmap$macrophage_sd = apply(for_heatmap[,names(for_heatmap)=="macrophage"], 1, sd)
for_heatmap$monocyte_mean = apply(for_heatmap[,names(for_heatmap)=="monocyte"], 1, mean)
for_heatmap$monocyte_sd = apply(for_heatmap[,names(for_heatmap)=="monocyte"], 1, sd)

to_plot = gather(for_heatmap[,17:22], "Group", "AUC", c(3,5)) # check this

ggplot(to_plot, aes(x=Group,y=AUC,group=Gene)) + geom_point(shape=17) + geom_line(alpha=0.1) + facet_wrap(~Cluster) + theme_thesis(10)
