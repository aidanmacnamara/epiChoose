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


# BRING IN THP-1/U937 -----------------------------------------------------

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


# PLOT RESULTS ------------------------------------------------------------

sampleDists <- dist(t(assay(rld_cell_lines)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(rld_cell_lines$condition, rld_cell_lines$rep, sep="_")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, show_colnames=FALSE)

pca_and_plot(rld_cell_lines, annot_1=rld_cell_lines$condition, annot_2=rld_cell_lines$cell_line)


# DIFFERENTIAL PMA VS. BASELINE -------------------------------------------

design(dds_cell_lines)
dds_cell_lines = DESeq(dds_cell_lines)
res_cell_lines = tbl_df(results(dds_cell_lines, contrast=c("condition","PMA","Baseline")))
res_cell_lines$gene = rownames(res_cell_lines)
res_cell_lines$symbol = gene_list_all$hgnc_symbol
arrange(res_cell_lines, padj)


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


# WHAT IS DRIVING PC1 DIFFERENCES? ----------------------------------------

y = t(assays(rld_all)[[1]])
dim(y)
pca_res <- prcomp(y, scale=TRUE, center=TRUE)

dat = data.frame(gene=gene_list_all$hgnc_symbol, loadings=pca_res$rotation[,1])
ranks = deframe(dat)
p = gmtPathways("tmp/c2.cp.v6.2.symbols.gmt")

res <- fgsea(pathways=p, stats=ranks, nperm=1000)
res_tidy = res %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()

print(ggplot(filter(res_tidy, padj<0.03), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.03)) + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score") + theme_thesis(10))

# remove the first pc and project again?
y_rev = pca_res$x[,-1] %*% t(pca_res$rotation[,-1])
dim(y_rev)
pca_and_plot(y_rev, annot_1=paste(rld_all$cell_type, rld_all$treatment, rld_all$time, sep="_"), annot_2=NA) # still orthogonal, easier to look gene-by-gene


# FC TABLE ----------------------------------------------------------------

# sanquin

design(dds_sanquin)
dds_sanquin$time_treatment = factor(paste(dds_sanquin$time, dds_sanquin$treatment, sep="_"))
design(dds_sanquin) = formula(~time_treatment+donor)
dds_sanquin = DESeq(dds_sanquin)

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

design(dds_cell_lines)
dds_cell_lines$cell_condition = factor(paste(dds_cell_lines$cell_line,dds_cell_lines$condition,sep="_"))
design(dds_cell_lines) = formula(~cell_condition)
design(dds_cell_lines)

dds_cell_lines = DESeq(dds_cell_lines)
trts_thp1 = unique(dds_cell_lines$cell_condition)[2:6] # pick out the treatment conditions
trts_u937 = unique(dds_cell_lines$cell_condition)[8:12]
res_cell_lines = data.frame(genes=gene_list_all$hgnc_symbol)

# add contrasts across thp-1 and u937
for(i in 1:length(trts_thp1)) {
  res = results(dds_cell_lines, contrast=c("cell_condition",trts_thp1[i],"THP-1_Baseline"))
  res_cell_lines = cbind(res_cell_lines, res$log2FoldChange, res$padj)
}

for(i in 1:length(trts_u937)) {
  res = results(dds_cell_lines, contrast=c("cell_condition",trts_u937[i],"U937_Baseline"))
  res_cell_lines = cbind(res_cell_lines, res$log2FoldChange, res$padj)
}

names(res_cell_lines)[-1] = paste(rep(c(trts_thp1,trts_u937), each=2), c("fc","p"), sep="_")
res_cell_lines = gather(res_cell_lines, "group","score", 2:dim(res_cell_lines)[2])
res_cell_lines$type = str_extract(res_cell_lines$group, "[a-z]+$")
res_cell_lines$group = str_replace(res_cell_lines$group, "_[a-z]+$", "")

res_all = rbind(res_cell_lines, res_sanquin) # fc and p values
# split fc and p values to 2 data frames
res_all_p = filter(res_all, type=="p") %>% dplyr::select(-type) %>% spread(group, score)
res_all_fc = filter(res_all, type=="fc") %>% dplyr::select(-type) %>% spread(group, score)
rownames(res_all_p) = res_all_p$genes
rownames(res_all_fc) = res_all_fc$genes
res_all_p = res_all_p[,-1]
res_all_fc = res_all_fc[,-1]
all(names(res_all_p)==names(res_all_fc))


# WHAT GENE SETS TO USE? --------------------------------------------------

# my_p = gmtPathways("tmp/c7.all.v6.2.symbols.gmt") # immune set
my_p = gmtPathways("tmp/c5.bp.v6.2.symbols.gmt") # biological processes
my_p_filter = my_p[grep("monocyte|macrophage", names(my_p), ignore.case=TRUE)]
my_p_filter = c(my_p_filter, my_p[grep("brain", names(my_p), ignore.case=TRUE)])
my_p_annot = data.frame(name=names(my_p_filter), tissue=c(rep("myeloid",20),rep("brain",14)), mean_fc=NA, sig_fcs=NA)

# my_p_filter = my_p_filter[1:20]
# my_p_annot = my_p_annot[1:20,]

lapply(my_p_filter, length)
my_p_filter_tbl = as.data.frame.matrix((table(stack(my_p_filter))))
my_p_filter_tbl = cbind(rownames(my_p_filter_tbl), my_p_filter_tbl)
rownames(my_p_filter_tbl) = NULL
upset(my_p_filter_tbl, order.by="freq", nsets=20, text.scale=1.5)

# CORRELATION BETWEEN CELL LINES AND PRIMARY ------------------------------

# heatmap annotation
s_annot = data.frame(Group=c(rep("SANQUIN",12),rep("GSK",10)), row.names=c(trts_sanquin,trts_thp1,trts_u937))

# store correlation between thp and sanquin
# cell_line_names = c(trts_thp1, trts_u937)
cell_line_names = c(trts_thp1, trts_u937, "5days_BG") # add control
sanquin_names = c("5days_BG","6days_Naive","5days_LPS","24hrs_BG","24hrs_Naive","24hrs_LPS","4hrs_BG","4hrs_Naive","4hrs_LPS","1hr_BG","1hr_Naive","1hr_LPS")
res_template = matrix(NA, nrow=length(cell_line_names), ncol=length(sanquin_names))
colnames(res_template) = sanquin_names
rownames(res_template) = cell_line_names

# store the results for each pathway
cor_dat = vector("list", length(my_p_filter))
names(cor_dat) = names(my_p_filter)

for(i in 1:length(cor_dat)) {
  
  res_copy = res_template
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]]) # get the gene indices for the pathway
  if(is_empty(g_ix)) next
  pathway_fc = res_all_fc[g_ix,] # get the fold changes
  my_p_annot$mean_fc[i] = mean(as.numeric(unlist(pathway_fc)), na.rm=TRUE)
  pathway_p = res_all_p[g_ix,]
  pathway_fc[pathway_p > 0.05] = NA # set non-significant (padj) to na
  my_p_annot$sig_fcs[i] = sum(pathway_p<0.05, na.rm=TRUE) / length(unlist(pathway_p))
  
  for(j in 1:length(cell_line_names)) {
    for(k in 1:length(sanquin_names)) {
      
      x = pathway_fc[,which(names(pathway_fc)==cell_line_names[j])]
      y = pathway_fc[,which(names(pathway_fc)==sanquin_names[k])]
      # plot(x,y)
      
      fisher_mat = matrix(
        c(
          sum(x>0 & y>0, na.rm=TRUE),
          sum(x<0 & y>0, na.rm=TRUE),
          sum(x>0 & y<0, na.rm=TRUE),
          sum(x<0 & y<0, na.rm=TRUE)
        ),
        nrow=2,
        dimnames=list(cell_line=c("up","down"),primary_cell=c("up","down"))
      )
      
      fisher_res = fisher.test(fisher_mat, alternative="greater")
      res_copy[j,k] = fisher_res$p.value
      
    }
  }
  
  cor_dat[[i]] = res_copy
  # pheatmap(pathway_fc, annotation_col=s_annot, main=names(my_p)[i], fontsize_row=5)
  # dat_tree = hclust(dist(t(dat_fc)))
  
}

ggplot(my_p_annot, aes(x=tissue, y=mean_fc)) + geom_boxplot() + theme_thesis(20) + xlab("Tissue") + ylab("Mean FC per Pathway")
ggplot(my_p_annot, aes(x=tissue, y=sig_fcs)) + geom_boxplot() + theme_thesis(15) + xlab("Tissue") + ylab("Proportion significant FCs per Pathway")

# make cor_dat long to look at correlation sets
cor_dat_long = data.frame()
my_labels = as.character(outer(cell_line_names, sanquin_names, paste, sep="_VS_"))
for(i in 1:length(cor_dat)) {
  print(i)
  cor_dat_long = rbind(
    cor_dat_long,
    data.frame(
      `P-value` = as.numeric(cor_dat[[i]]),
      Label = my_labels,
      Pathway = names(cor_dat)[i], check.names=FALSE
    )
  )
}
cor_dat_long$tissue = my_p_annot$tissue[match(cor_dat_long$Pathway, my_p_annot$name)]
cor_dat_long = tbl_df(cor_dat_long)

# WHAT ARE THE HITS? ------------------------------------------------------

# the overall distribution of fisher exact test p-values
cor_dat_long %>% ggplot(aes(x=`P-value`)) + geom_histogram(binwidth=0.01, color="black") + theme_thesis()

# are related primary cells more highly correlated across myeloid pathways compared to brain pathways?
cor_dat_long %>% filter(Label=="5days_BG_VS_6days_Naive") %>% ggplot(aes(x=tissue, y=`P-value`)) + geom_boxplot() + theme_thesis()

cutoff = 0.05
cor_dat_long %>% filter(`P-value` < cutoff) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N))
cor_dat_long %>% filter(`P-value` < cutoff) %>% group_by(Label) %>% summarise(N=n()) %>% arrange(desc(N))
cor_dat_long %>% filter(`P-value` < cutoff, grepl("U937", Label), tissue=="myeloid")

select_ps = cor_dat_long %>% filter(`P-value` < cutoff, grepl("THP-1|U937",Label)) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N)) %>% select(Pathway) %>% unlist %>% as.character()

for(j in 1:length(select_ps)) {
  
  i = which(names(my_p_filter)==select_ps[j])
  i = which(names(my_p_filter)=="GO_MONOCYTE_CHEMOTAXIS")
  cor_dat[[i]]
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]])
  dat_fc = res_all_fc[g_ix,]
  pheatmap(dat_fc, main=names(my_p_filter)[i], fontsize_row=5)
  
}

gsea_in = res_all_fc$`5days_LPS`
names(gsea_in) = rownames(res_all_fc)
res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()


# FUNCTIONS ---------------------------------------------------------------

pca_and_plot <- function(rld, annot_1, annot_2) {
  
  if(class(rld)=="matrix") {
    y = rld
  } else {
    y = t(assays(rld)[[1]])
  }
  
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


