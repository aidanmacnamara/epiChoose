# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanqin primary data

load("data/total_reg.RData")
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
    dat_add = total_reg$H3K27ac$res[row_ix,which(roi_reg$feature_type_name=="Enhancer")]
  }
  if(names(dds_list)[i]=="atac") {
    dat_add = total_reg$ATAC$res[row_ix,which(roi_reg$feature_type_name=="Open chromatin")]
  }
  
  dat_add[is.na(dat_add)] = 0
  dat_add[dat_add==0] = 1
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