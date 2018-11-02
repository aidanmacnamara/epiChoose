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

png(filename="c:/Users/am673712/Dropbox/OTAR020/images/sanquin_image_2.png", width=1000, height=600)
pca_and_plot(rld_sanquin, annot_1=paste(rld_sanquin$treatment, rld_sanquin$time, sep="_"), annot_2=rld_sanquin$donor)
dev.off()


# TOP N HITS --------------------------------------------------------------

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

# PICK SAMPLES ------------------------------------------------------------

row_ix = which(grepl("thp-1", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE))

col_data = dat_all$tss$H3K27ac$annot[row_ix,]
col_data_filt = data.frame(dplyr::select(col_data, Label))
col_data_filt$rep = str_extract(col_data_filt$Label, "BR[12]")
col_data_filt$condition = str_extract(col_data_filt$Label, "[[:alnum:]+]+$")


# RUN DDS -----------------------------------------------------------------

dat_add = dat_all$tss$H3K27ac$res[row_ix,]
# dat_add = total_reg$H3K27ac$res[row_ix,which(roi_reg$feature_type_name=="Enhancer")]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_thp = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~rep+condition)
rld_thp = vst(dds_thp, blind=FALSE)


# PLOT RESULTS ------------------------------------------------------------

sampleDists <- dist(t(assay(rld_thp)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_thp$condition, rld_thp$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_thp, annot_1=rld_thp$condition, annot_2=NA)


# DIFFERENTIAL PMA VS. BASELINE -------------------------------------------

design(dds_thp)
dds_thp = DESeq(dds_thp)
res_thp = tbl_df(results(dds_thp, contrast=c("condition","PMA","Baseline")))
res_thp$gene = rownames(res_thp)
res_thp$symbol = gene_list_all$hgnc_symbol
arrange(res_thp, padj)


# COMPARE THP1 AND SANQUIN TOGETHER ---------------------------------------

# PICK SAMPLES ------------------------------------------------------------

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

pca_and_plot(rld_all, annot_1=paste(rld_all$treatment, rld_all$time, sep="_"), annot_2=rld_all$donor)


# WHAT IS DRIVING PC1 DIFFERENCES? ----------------------------------------

dat = data.frame(gene=gene_list_all$hgnc_symbol, loadings=pca_res$rotation[,1])
ranks = deframe(dat)
p = gmtPathways("../sp_140_follow_up/c2.cp.v6.2.symbols.gmt")

res <- fgsea(pathways=p, stats=ranks, nperm=1000)
res_tidy = res %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()

print(ggplot(filter(res_tidy, padj<0.03), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.03)) + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score") + theme_thesis(10))

# remove the first pc and project again?
y_rev = pca_res$x[,-1] %*% t(pca_res$rotation[,-1])
pca_and_plot(y_rev) # still orthogonal, easier to look gene-by-gene


# SEARCH BY GENE SETS -----------------------------------------------------

my_p = gmtPathways("tmp/c2.cp.reactome.v6.2.symbols.gmt")


# FC TABLE ----------------------------------------------------------------

# sanquin

design(dds_sanquin)
trts_sanquin = as.character(unique(dds_sanquin$time_treatment)[-1])
res_sanquin = data.frame(genes=gene_list_all$hgnc_symbol)
for(i in 1:length(trts_sanquin)) {
  res = results(dds_sanquin, contrast=c("time_treatment",trts_sanquin[i],"1hr_Attached"))
  res_sanquin = cbind(res_sanquin, res$log2FoldChange, res$padj)
}
names(res_sanquin)[-1] = paste(rep(trts_sanquin, each=2), c("fc","p"), sep="_")
res_sanquin = gather(res_sanquin, "group","score", 2:dim(res_sanquin)[2])
res_sanquin$type = str_extract(res_sanquin$group, "[a-z]+$")
res_sanquin$group = str_replace(res_sanquin$group, "_[a-z]+$", "")

# gsk

design(dds_thp)
trts_thp = as.character(unique(dds_thp$condition)[-1])
res_thp = data.frame(genes=gene_list_all$hgnc_symbol)
for(i in 1:length(trts_thp)) {
  res = results(dds_thp, contrast=c("condition",trts_thp[i],"Baseline"))
  res_thp = cbind(res_thp, res$log2FoldChange, res$padj)
}
names(res_thp)[-1] = paste(rep(trts_thp, each=2), c("fc","p"), sep="_")
res_thp = gather(res_thp, "group","score", 2:dim(res_thp)[2])
res_thp$type = str_extract(res_thp$group, "[a-z]+$")
res_thp$group = str_replace(res_thp$group, "_[a-z]+$", "")

res_all = rbind(res_thp, res_sanquin) # fc and p values
# split fc and p values to 2 data frames
res_all_p = filter(res_all, type=="p") %>% dplyr::select(-type) %>% spread(group, score)
res_all_fc = filter(res_all, type=="fc") %>% dplyr::select(-type) %>% spread(group, score)
rownames(res_all_p) = res_all_p$genes
rownames(res_all_fc) = res_all_fc$genes
res_all_p = res_all_p[,-1]
res_all_fc = res_all_fc[,-1]
all(names(res_all_p)==names(res_all_fc))

# heatmap annotation
s_annot = data.frame(Group=c(rep("GSK",5),rep("SANQUIN",12)), row.names=c(trts_thp,trts_sanquin))

# store correlation between thp and sanquin
cor_dat = vector("list", length(my_p))
names(cor_dat) = names(my_p)
thp_names = c("LPS","PMA","VD3","VD3+LPS","PMA+LPS")
sanquin_names = c("5days_BG","6days_Naive","5days_LPS")
res_template = matrix(NA, nrow=length(thp_names), ncol=length(sanquin_names))
colnames(res_template) = sanquin_names
rownames(res_template) = thp_names

i=7
for(i in 1:10) {
  
  res_copy = res_template
  g_ix = which(rownames(res_all_fc) %in% my_p[[i]])
  if(is_empty(g_ix)) next
  dat_fc = res_all_fc[g_ix,]
  for(j in 1:5) {
    for(k in 1:3) {
      y = cor.test(
        dat_fc[,which(names(dat_fc)==thp_names[j])],
        dat_fc[,which(names(dat_fc)==sanquin_names[k])]
      )
      # if(y$p.value<0.05) {
        res_copy[j,k] = y$estimate
        
      # }
    }
  }
  
  cor_dat[[i]] = res_copy
  # dat_tree = hclust(dist(t(dat_fc)))
  # pheatmap(dat_fc, annotation_col=s_annot, main=names(my_p)[i], fontsize_row=5)
  # mix_dat$mixing[i] = all(c(1,2) %in% cutree(dat_tree, k=2)[s_annot$Group=="GSK"])
  
}

table(mix_dat$mixing)
head(which(mix_dat$mixing))

gsea_in = res_all_fc$`5days_LPS`
names(gsea_in) = rownames(res_all_fc)
res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()


# FUNCTIONS ---------------------------------------------------------------

pca_and_plot <- function(rld, annot_1, annot_2) {
  
  y = t(assays(rld)[[1]])
  dim(y)
  y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(y)
  
  pca_res <- prcomp(y, scale=TRUE, center=TRUE)
  pca_res_summary = summary(pca_res)
  yy = data.frame(pca_res$x[,1:2])
  names(yy) = c("x","y")
  yy$annot_1 = annot_1
  yy$annot_2 = annot_2
  
  if(is.na(annot_2)) {
    my_plot = ggplot(yy, aes(x=x, y=y)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")
  } else {
    my_plot = ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")
  }
  return(my_plot)
  
}


