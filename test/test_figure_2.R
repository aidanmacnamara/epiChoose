# comparison of monocytes/macrophages with u937/thp-1

# 1. establish the amount of paired samples in blueprint
# 2. groups: monocytes, macrophages, u937, thp-1
# 3. how to define dynamic regions?
# 4. figure 4d for dynamic regions i.e. regions that differ between monocytes and macrophages and what they look like in u937/thp-1
# 5. tf inference from deregulated sites (look at luz's review for methods)

load("data/total_data.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")


# ESTABLISH DATA MATRICES -------------------------------------------------

# enhancer/promoter/atac matrices
dds_list = vector("list", 3)
names(dds_list) = c("promoter","enhancer","atac")
rld_list = vector("list", 3)
names(rld_list) = c("promoter","enhancer","atac")

# pick samples
row_ix = which(grepl("thp-1", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("monocyte", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("macrophage", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE) | grepl("u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE))

col_data = dat_all$tss$H3K27ac$annot[row_ix,]
col_data_filt = data.frame(select(col_data, Label))
col_data_filt$cell = tolower(str_extract(col_data_filt$Label, "[Mm]onocyte|macrophage|THP-1|U937"))
col_data_filt$donor = str_extract(col_data_filt$Label, "^[CS][[:alnum:]]+")
col_data_filt$condition = tolower(str_extract(col_data_filt$Label, "inflammatory|classical|alternatively activated|BR.*$"))
col_data_filt$condition = str_replace(col_data_filt$condition, "br[12]_", "")
col_data_filt$condition[is.na(col_data_filt$condition)] = "primary"
rownames(col_data_filt) = rownames(dat_all$tss$H3K27ac$res)[row_ix]

roi_reg_enhancer = roi_reg[which(roi_reg$feature_type_name=="Open chromatin")]

for(i in 1:length(dds_list)) {
  
  if(names(dds_list)[i]=="promoter") {
    dat_add = dat_all$tss$H3K27ac$res[row_ix,]
  }  
  if(names(dds_list)[i]=="enhancer") {
    dat_add = total_data$H3K27ac$res[row_ix,which(roi_reg$feature_type_name=="Enhancer")]
  }
  if(names(dds_list)[i]=="atac") {
    dat_add = total_data$ATAC$res[row_ix,which(roi_reg$feature_type_name=="Open chromatin")]
  }
  
  dat_add[is.na(dat_add)] = 0
  dat_add = apply(dat_add, 2, as.integer)
  dim(dat_add)
  dds = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~condition)
  dds_list[[i]] = dds
  rld_list[[i]] = rlog(dds, blind=FALSE)  
  
}

for(i in 1:2) {
  
  sampleDists <- dist(t(assay(rld_list[[i]])))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld_list[[i]]$cell, rld_list[[i]]$condition, sep="_")
  # colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  png(filename=paste0("c:/Users/am673712/Dropbox/OTAR020/images/image_1_", i, ".png"), width=1192, height=652)
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)
  dev.off()
  
  # pca
  
  y = t(assays(rld_list[[i]])[[1]])
  dim(y)
  y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(y)
  
  pca_res <- prcomp(y, scale=TRUE, center=TRUE)
  pca_res_summary = summary(pca_res)
  yy = data.frame(pca_res$x[,1:2])
  names(yy) = c("x","y")
  yy$annot_1 = paste(colData(rld_list[[i]])$cell, colData(rld_list[[i]])$condition, sep="_")
  png(filename=paste0("c:/Users/am673712/Dropbox/OTAR020/images/image_2_", i, ".png"), width=1430, height=788)
  print(ggplot(yy, aes(x=x, y=y)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5))
  dev.off()
  
  design(dds_list[[i]])
  dds_list[[i]] = DESeq(dds_list[[i]], test="LRT", reduced=~1)
  res = tbl_df(results(dds_list[[i]]))
  res$gene = rownames(results(dds_list[[i]]))
  res
  
  for(j in 1:3) {
    y = plotCounts(dds_list[[i]], gene=head(order(res$padj),10)[j], intgroup="condition", returnData=TRUE)
    png(filename=paste0("c:/Users/am673712/Dropbox/OTAR020/images/image_3_", i, "_", j, ".png"), width=1106, height=569)
    print(ggplot(y, aes(x=condition, y=count)) + geom_jitter(shape=17, size=3) + theme_thesis(20) + xlab("") + ylab("Count") + ggtitle(res$gene[head(order(res$padj),10)[j]]))
    dev.off()
  }
  
  # take the top 1000 hits
  genes_top = head(res$gene[order(res$padj)],1e3)
  
  for_heatmap = assays(rld_list[[i]])[[1]]
  colnames(for_heatmap) = paste(col_data_filt$cell, col_data_filt$condition, sep="_")
  for_heatmap = as.data.frame(for_heatmap[match(genes_top,rownames(dds_list[[i]])),])
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
  
  png(filename=paste0("c:/Users/am673712/Dropbox/OTAR020/images/image_4_", i, ".png"), width=794, height=1287)
  pheatmap(for_heatmap, fontsize_row=8,
           annotation_row = data.frame(
             K_Means = factor(fit$cluster),
             row.names=rownames(for_heatmap)
           ), show_rownames=FALSE
  )
  dev.off()
  
  for_heatmap$Gene = rownames(for_heatmap)
  for_heatmap$Cluster = factor(fit$cluster)
  head(for_heatmap)
  for_heatmap = gather(for_heatmap, "Group", "AUC", 1:17)
  for_heatmap$Group = factor(for_heatmap$Group)
  png(filename=paste0("c:/Users/am673712/Dropbox/OTAR020/images/image_5_", i, ".png"), width=1089, height=748)
  ggplot(for_heatmap, aes(x=Group,y=AUC,group=Gene)) + geom_point(shape=17) + geom_line(alpha=0.1) + facet_wrap(~Cluster) + theme_thesis(15)
  dev.off()
  
}


# MONOCYTE/MACROPHAGE RELEVANT GENES --------------------------------------

# pick the gene set from saeed/novakovic
p_genes = read_excel("tmp/1-s2.0-S0092867416313162-mmc2.xlsx", sheet="Table S2B. Genes H3K27ac prom", skip=1)
dg_genes = p_genes$Promoter_H3K27ac_Differentiation_gain; dg_genes = dg_genes[!is.na(dg_genes)]
dl_genes = p_genes$Promoter_H3K27ac_Differentiation_loss; dl_genes = dl_genes[!is.na(dl_genes)]

for_heatmap = assays(rld_list[[1]])[[1]]
colnames(for_heatmap) = col_data_filt$cell
match_ix = match(c(dg_genes,dl_genes), gene_list_all$hgnc_symbol); match_ix = match_ix[!is.na(match_ix)]

for_heatmap = as.data.frame(for_heatmap[match_ix,grepl("monocyte|macrophage", colnames(for_heatmap))])
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


# COUNT AUC COMPARISON ----------------------------------------------------

# regions

load("data/dat_all.RData")
tss_window = 2e3
tss_regions = gene_list_all
start(tss_regions) = tss_regions$transcription_start_site - tss_window
end(tss_regions) = tss_regions$transcription_start_site + tss_window

# counts

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")
bamfiles = filter(data_gsk, grepl("thp", Label, ignore.case=TRUE), Mark=="H3K27ac") %>% dplyr::select(Bam) %>% unlist()

col_data = data_gsk[match(str_extract(bamfiles, "[[:alnum:]\\._-]+$"), str_extract(data_gsk$Bam, "[[:alnum:]\\._-]+$")),] %>% dplyr::select(Cell, Mark, Rep, Stimulus)
col_data = as.data.frame(col_data)

bamfiles <- BamFileList(bamfiles, yieldSize=2000000)
lapply(bamfiles, seqinfo)
register(MulticoreParam())
se <- summarizeOverlaps(features=tss_regions, reads=bamfiles, mode="Union", ignore.strand=TRUE)

rownames(col_data) = rownames(colData(se))
colData(se) <- DataFrame(col_data)

dds = DESeqDataSet(se, design=~Stimulus)
dds = DESeq(dds)
rld = rlog(dds, blind=FALSE) 

res = results(dds, contrast=c("Stimulus","PMA","Baseline"))
res$hgnc_symbol = gene_list_all$hgnc_symbol
res_count = tbl_df(res)
save(res_count, file="tmp/res_count.RData")

# auc with limma

load("tmp/fc_res_long.RData")
res_auc = filter(fc_res_long, trmt=="THP_1_PMA")

to_plot_comparison = data.frame(
  genes = gene_list_all$hgnc_symbol,
  `LogFC Count` = res_count$log2FoldChange[match(gene_list_all$hgnc_symbol, res_count$hgnc_symbol)],
  `LogFC AUC` = res_auc$fc[match(gene_list_all$hgnc_symbol, res_auc$gene)],
  Group = NA, check.names=FALSE
)

cutoff = 1
to_plot_comparison$Group[to_plot_comparison$genes %in% intersect(res_auc$gene[abs(res_auc$fc)>=cutoff & res_auc$p<=0.05], res_count$hgnc_symbol[abs(res_count$log2FoldChange)>=cutoff & res_count$padj<=0.05])] = "Both"
to_plot_comparison$Group[to_plot_comparison$genes %in% setdiff(res_auc$gene[abs(res_auc$fc)>=cutoff & res_auc$p<=0.05], res_count$hgnc_symbol[abs(res_count$log2FoldChange)>=cutoff & res_count$padj<=0.05])] = "AUC"
to_plot_comparison$Group[to_plot_comparison$genes %in% setdiff(res_count$hgnc_symbol[abs(res_count$log2FoldChange)>=cutoff & res_count$padj<=0.05], res_auc$gene[abs(res_auc$fc)>=cutoff & res_auc$p<=0.05])] = "Count"
to_plot_comparison$Group[is.na(to_plot_comparison$Group)] = "Not Significant"
to_plot_comparison$size = ifelse(to_plot_comparison$Group=="Not Significant", 0.1, 1)

ggplot(to_plot_comparison, aes(x=`LogFC Count`, y=`LogFC AUC`, color=factor(Group))) + geom_point(aes(size=factor(size)), shape=17, alpha=0.5) + theme_thesis() + scale_size_discrete(guide=FALSE)
save(to_plot_comparison, file="tmp/to_plot_comparison.RData")


