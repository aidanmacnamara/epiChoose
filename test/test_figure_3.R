
# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanquin primary data

# LOAD --------------------------------------------------------------------

require(tidyverse)
require(RColorBrewer)
require(ggrepel)
require(RcisTarget)

loadload("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")

load("tmp/dds_sanquin.RData")
load("tmp/dds_cell_lines.RData")

load("tmp/rld_sanquin.RData")
load("tmp/rld_cell_lines.RData")


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


# WHAT GENE SETS TO USE? --------------------------------------------------

# my_p = gmtPathways("tmp/c7.all.v6.2.symbols.gmt") # immune set
# my_p = gmtPathways("tmp/c2.cp.reactome.v6.2.symbols.gmt") # biological processes
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


# FC TABLE ----------------------------------------------------------------

# sanquin

design(dds_sanquin)
dds_sanquin$time_treatment = factor(paste(dds_sanquin$time, dds_sanquin$treatment, sep="_"))
design(dds_sanquin) = formula(~time_treatment+donor)
dds_sanquin = DESeq(dds_sanquin)

trts_sanquin = as.character(unique(dds_sanquin$time_treatment)[-1])
res_sanquin = data.frame(genes=gene_list_all$hgnc_symbol)
res_sanquin_plot = data.frame(trmt=NULL, fc=NULL, p=NULL)

for(i in 1:length(trts_sanquin)) {
  res = results(dds_sanquin, contrast=c("time_treatment",trts_sanquin[i],"1hr_Attached"))
  res_sanquin = cbind(res_sanquin, res$log2FoldChange, res$padj)
  res_sanquin_plot = rbind(res_sanquin_plot, data.frame(trmt=trts_sanquin[i], fc=res$log2FoldChange, p=res$padj))
}

res_sanquin_plot$gene = rep(gene_list_all$hgnc_symbol, each=12)
names(res_sanquin)[-1] = paste(rep(trts_sanquin, each=2), c("fc","p"), sep="_")

# try to replicate figure 1b from novakovic here
# filter for dynamic regions across all conditions
dyn_genes = filter(res_sanquin_plot, p <= 0.05, abs(fc) >= 3) %>% dplyr::select(gene) %>% unlist() %>% as.character()
pca_and_plot(rld_sanquin[which(gene_list_all$hgnc_symbol %in% dyn_genes)], annot_1=paste(rld_sanquin$treatment, rld_sanquin$time, sep="_"), annot_2=rld_sanquin$donor)

res_sanquin = gather(res_sanquin, "group","score", 2:dim(res_sanquin)[2])
res_sanquin$type = str_extract(res_sanquin$group, "[a-z]+$")
res_sanquin$group = str_replace(res_sanquin$group, "_[a-z]+$", "")

# gsk

design(dds_cell_lines)
dds_cell_lines$cell_condition = factor(paste(dds_cell_lines$cell_line,dds_cell_lines$condition,sep="_"))
design(dds_cell_lines) = formula(~cell_condition)
design(dds_cell_lines)

dds_cell_lines = DESeq(dds_cell_lines)
trts_thp1 = as.character(unique(dds_cell_lines$cell_condition)[2:6]) # pick out the treatment conditions
trts_u937 = as.character(unique(dds_cell_lines$cell_condition)[8:12])
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


# CONDITION ENRICHMENT ----------------------------------------------------

# what is the enrichment for a condition, e.g. thp1 + pma

conds = c(trts_thp1,trts_u937)
conds_list = vector("list", length(conds))
names(conds_list) = conds

for(i in 1:length(conds)) {
  
  my_cond = res_cell_lines[,c(1, which(str_replace(names(res_cell_lines), "_[pfc]+$", "") %in% conds[i]))]
  my_cond = my_cond[my_cond[,3] < 0.05,]
  my_cond = my_cond[order(-my_cond[,2]),]
  my_cond = my_cond[!is.na(my_cond$genes),]
  my_cond %>% head()
  
  # gsea_in = my_cond[,2]
  # names(gsea_in) = my_cond$genes
  # res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
  # res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
  # res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()
  
  conds_list[[i]] = my_cond
  
}

gene_lists = lapply(conds_list, function(x) as.character(x[x[,2]>2, 'genes']))

# import motif-gene rankings
# how is this ranking performed?
# how much do the results change varying the tss window?
motif_rankings <- importRankings("tmp/hg19-tss-centered-5kb-7species.mc9nr.feather")

# import motif-tf annotations
data(motifAnnotations_hgnc)

# run enrichment
# question: what motifs are enriched in each cluster?
motif_enrichment <- cisTarget(gene_lists, motif_rankings, motifAnnot=motifAnnotations_hgnc)

for(j in 1:length(gene_lists)) { # visualise
  
  # get the significant motifs
  # can also apply motif-tf  mapping threshold here
  sig_motifs = filter(motif_enrichment, geneSet==conds[j], NES>=10) %>% select(motif) %>% unlist()
  
  # get all interactions between significant tfs and regulated genes
  inc_mat <- getSignificantGenes(gene_lists[[j]], motif_rankings, signifRankingNames=sig_motifs, plotCurve=FALSE, maxRank=5000-20, genesFormat="incidMatrix", method="aprox")$incidMatrix
  
  # construct network  
  edges <- melt(inc_mat)
  edges <- edges[which(edges[,3]==1),1:2]
  colnames(edges) <- c("from","to")
  motifs <- unique(as.character(edges[,1]))
  genes <- unique(as.character(edges[,2]))
  nodes <- data.frame(id=c(motifs, genes), label=c(motifs, genes), title=c(motifs, genes), shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))), color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
  
  # visualise
  visNetwork(nodes, edges) %>% visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>% visPhysics(enabled=FALSE)
  
}


# \ CONDITION ENRICHMENT --------------------------------------------------

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


# CORRELATION BETWEEN CELL LINES AND PRIMARY ------------------------------

# store correlation between thp and sanquin
trts_cell_lines = c(trts_thp1, trts_u937, "5days_BG") # add 5 days bg as control
res_template = matrix(NA, nrow=length(trts_cell_lines), ncol=length(trts_sanquin))
colnames(res_template) = trts_sanquin
rownames(res_template) = trts_cell_lines

# store the results for each pathway
out_comp = vector("list", length(my_p_filter))
names(out_comp) = names(my_p_filter)

out_activity = matrix(NA, nrow=length(my_p_filter), ncol=length(c(trts_cell_lines, trts_sanquin)))
colnames(out_activity) = c(trts_cell_lines, trts_sanquin)
rownames(out_activity) = names(my_p_filter)

for(i in 1:length(out_comp)) {
  
  res_copy = res_template
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]]) # get the gene indices for the pathway
  if(is_empty(g_ix)) next
  
  pathway_fc = res_all_fc[g_ix,] # get the fold changes
  pathway_p = res_all_p[g_ix,] # get the p-values
  my_p_annot$mean_fc[i] = mean(as.numeric(unlist(pathway_fc)), na.rm=TRUE) # what is the mean fc for the pathway?
  my_p_annot$sig_fcs[i] = sum(pathway_p<0.05, na.rm=TRUE) / length(unlist(pathway_p)) # how many sig fcs for the pathway?
  
  # pathway_fc[pathway_p > 0.05] = NA # set non-significant (padj) fcs to na
  
  fcs_per_condition = apply(pathway_fc, 2, function(x) sum(x < 0.05) / dim(pathway_fc)[1])
  out_activity[i,] = fcs_per_condition[match(colnames(out_activity),names(fcs_per_condition))]
  
  for(j in 1:length(trts_cell_lines)) {
    for(k in 1:length(trts_sanquin)) {
      
      x = pathway_fc[,which(names(pathway_fc)==trts_cell_lines[j])]
      y = pathway_fc[,which(names(pathway_fc)==trts_sanquin[k])]
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
  
  out_comp[[i]] = res_copy
  
}

ggplot(my_p_annot, aes(x=tissue, y=mean_fc)) + geom_boxplot() + theme_thesis(20) + xlab("Tissue") + ylab("Mean FC per Pathway")
ggplot(my_p_annot, aes(x=tissue, y=sig_fcs)) + geom_boxplot() + theme_thesis(15) + xlab("Tissue") + ylab("Proportion significant FCs per Pathway")

# make out_comp long to look at correlation sets
out_comp_long = data.frame()
my_labels = as.character(outer(trts_cell_lines, trts_sanquin, paste, sep="_VS_"))
for(i in 1:length(out_comp)) {
  print(i)
  out_comp_long = rbind(
    out_comp_long,
    data.frame(
      `P-value` = as.numeric(out_comp[[i]]),
      Label = my_labels,
      Pathway = names(out_comp)[i], check.names=FALSE
    )
  )
}
out_comp_long$tissue = my_p_annot$tissue[match(out_comp_long$Pathway, my_p_annot$name)]
out_comp_long = tbl_df(out_comp_long)


# CONDITIONS VS PATHWAYS --------------------------------------------------

pheatmap(out_activity)


# WHAT ARE THE HITS? ------------------------------------------------------

# the overall distribution of fisher exact test p-values
out_comp_long %>% ggplot(aes(x=`P-value`)) + geom_histogram(binwidth=0.01, color="black") + theme_thesis()

# are related primary cells more highly correlated across myeloid pathways compared to brain pathways?
out_comp_long %>% filter(Label=="5days_BG_VS_6days_Naive") %>% ggplot(aes(x=tissue, y=`P-value`)) + geom_boxplot() + theme_thesis()

cutoff = 0.05
out_comp_long %>% filter(`P-value` < cutoff) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N))
out_comp_long %>% filter(`P-value` < cutoff) %>% group_by(Label) %>% summarise(N=n()) %>% arrange(desc(N))
out_comp_long %>% filter(`P-value` < cutoff, grepl("U937", Label), tissue=="myeloid")

select_ps = out_comp_long %>% filter(`P-value` < cutoff, grepl("THP-1|U937",Label)) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N)) %>% select(Pathway) %>% unlist %>% as.character()

for(j in 1:length(select_ps)) {
  
  # i = which(names(my_p_filter)==select_ps[j])
  i = which(names(my_p_filter)=="GO_MACROPHAGE_CHEMOTAXIS")
  out_comp[[i]]
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]])
  dat_fc = res_all_fc[g_ix,]
  pheatmap(dat_fc, main=names(my_p_filter)[i], fontsize_row=5)
  
}

gsea_in = res_all_fc$`5days_LPS`
names(gsea_in) = rownames(res_all_fc)
res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()


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


# DIFFERENTIAL PMA VS. BASELINE -------------------------------------------

design(dds_cell_lines)
dds_cell_lines = DESeq(dds_cell_lines)
res_cell_lines = tbl_df(results(dds_cell_lines, contrast=c("condition","PMA","Baseline")))
res_cell_lines$gene = rownames(res_cell_lines)
res_cell_lines$symbol = gene_list_all$hgnc_symbol
arrange(res_cell_lines, padj)


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

