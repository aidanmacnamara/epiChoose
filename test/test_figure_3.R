
# look at specific comparisons

# LOAD --------------------------------------------------------------------

require(tidyverse)
require(RColorBrewer)
require(ggrepel)
require(RcisTarget)
require(DESeq2)
require(d3heatmap)
require(heatmaply)
require(gplots)
require(reshape2)
require(UpSetR)
require(eulerr)

load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")

load("tmp/dds_cell_lines.RData")
load("tmp/dds_saeed.RData")
load("tmp/dds_novakovic.RData")

load("tmp/rld_cell_lines.RData")
load("tmp/rld_saeed.RData")
load("tmp/rld_novakovic.RData")

load("tmp/rld_all.RData")


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


# GET FOLD CHANGES --------------------------------------------------------

novakovic = as.character(unique(dds_novakovic$time_treatment)[-1])
novakovic = novakovic[!grepl("^1hr", novakovic)]

trts = list(
  novakovic = novakovic,
  saeed = as.character(unique(dds_saeed$time_treatment))[1:3],
  cell_lines = c(
    as.character(unique(dds_cell_lines$cell_condition)[2:6]), # pick out the treatment conditions
    as.character(unique(dds_cell_lines$cell_condition)[8:12])
  )
)

fc_res = data.frame(genes=gene_list_all$hgnc_symbol)
fc_res_long = data.frame(trmt=NULL, fc=NULL, p=NULL, group=NULL)

for(i in 1:length(trts)) {
  for(j in 1:length(trts[[i]])) {
    
    if(i==1) {
      if(grepl("LPS", trts[[i]][j])) {
        res = results(dds_novakovic, contrast=c("time_treatment",trts[[i]][j],"1hr_LPS"))
      }
      if(grepl("BG", trts[[i]][j])) {
        res = results(dds_novakovic, contrast=c("time_treatment",trts[[i]][j],"1hr_BG"))
      }
      if(grepl("Naive", trts[[i]][j])) {
        res = results(dds_novakovic, contrast=c("time_treatment",trts[[i]][j],"1hr_Naive"))
      }
    }
    if(i==2) res = results(dds_saeed, contrast=c("time_treatment",trts[[i]][j],"0days_Untreated"))
    if(i==3) {
      if(grepl("THP-1", trts[[i]][j])) {
        res = results(dds_cell_lines, contrast=c("cell_condition",trts[[i]][j],"THP-1_Baseline"))
      } else {
        res = results(dds_cell_lines, contrast=c("cell_condition",trts[[i]][j],"U937_Baseline"))
      }
    }
    
    fc_res = cbind(fc_res, res$log2FoldChange, res$padj)
    fc_res_long = rbind(fc_res_long, data.frame(trmt=trts[[i]][j], fc=res$log2FoldChange, p=res$padj, group=names(trts)[i]))
  }
}

fc_res_long$gene = rep(gene_list_all$hgnc_symbol, length(unlist(trts)))
# names(fc_res)[-1] = paste(c(rep(names(trts)[1],length(trts[[1]])*2), rep(names(trts)[2],length(trts[[2]])*2), rep(names(trts)[3],length(trts[[3]])*2)), paste(rep(unlist(trts), each=2), c("fc","p"), sep="_"), sep="_")
names(fc_res)[-1] = paste(rep(unlist(trts), each=2), c("fc","p"), sep="_")


# DYNAMIC REGIONS ---------------------------------------------------------

# define first across the primary data

fc_res_long$trmt = as.character(fc_res_long$trmt)

# fc_res_long_filt = filter(fc_res_long, p <= 0.05, abs(fc) >= 4)
fc_res_long_filt = filter(fc_res_long, p <= 0.05, abs(fc) >= 2, !grepl("BG|LPS|glucan|U937|THP", trmt, ignore.case=TRUE)) # filter out stimulatory conditions
fc_res_long_filt = arrange(fc_res_long_filt, desc(abs(fc)))

sort(table(fc_res_long_filt$trmt))
tail(sort(table(fc_res_long_filt$gene)), 20)

plot(euler(
  list(
    novakovic_4_hrs = unique(fc_res_long_filt[fc_res_long_filt$trmt=="4hrs_Naive", 'gene']),
    novakovic_24_hrs = unique(fc_res_long_filt[fc_res_long_filt$trmt=="24hrs_Naive", 'gene']),
    novakovic_6_days = unique(fc_res_long_filt[fc_res_long_filt$trmt=="6days_Naive", 'gene']),
    saeed_6_days = unique(fc_res_long_filt[fc_res_long_filt$trmt=="6days_Untreated", 'gene'])
  )
), quantities=TRUE)

dyn_genes = unique(fc_res_long_filt$gene) # primary dymamic h3k27ac signals

col_ix = match(dyn_genes, gene_list_all$hgnc_symbol)

y = assays(rld_all)[[1]]
dim(y)
rownames(y) = gene_list_all$hgnc_symbol
y = y[col_ix,]
dim(y)
y = y[apply(y, 1, function(x) sd(x)>1),]
dim(y)

y = y[,!grepl("BG|LPS|glucan", colnames(y), ignore.case=TRUE)]
dim(y)

d3heatmap(t(y), cexCol=0.8)
heatmaply(t(y), col=bluered(75))

filter(tbl_df(fc_res), genes %in% rownames(y)) %>% write_csv("out.csv")

# correlation matrix

yy = cor(y)
heatmaply(yy, col=bluered(75))

row_ix = grep("SANQUIN|N000", rownames(yy))
yyy = yy[row_ix, -row_ix]


# CONDITION ENRICHMENT ----------------------------------------------------

# what is the enrichment for a condition, e.g. thp1 + pma
conds = unlist(trts)
conds = conds[!grepl("BG|LPS|glucan", conds, ignore.case=TRUE)]

conds_list = vector("list", length(conds))
names(conds_list) = conds

tfs_list =  vector("list", length(conds))
names(tfs_list) = conds

for(i in 1:length(conds)) {
  
  my_cond = fc_res[,c(1, which(str_replace(names(fc_res), "_[pfc]+$", "") %in% conds[i]))]
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

gene_lists = lapply(conds_list, function(x) as.character(x[x[,2]>2, 'genes'])) # take everything with a fc > 2
lapply(gene_lists, length)

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
  
  print(j)
  
  # get the significant motifs
  # can also apply motif-tf  mapping threshold here
  sig_motifs = filter(motif_enrichment, geneSet==conds[j], NES>=4)
  
  if(dim(sig_motifs)[1]==0) next
  
  # get all interactions between significant tfs and regulated genes
  # this is also identifying the significant regulated genes
  inc_mat = getSignificantGenes(gene_lists[[j]], motif_rankings, signifRankingNames=sig_motifs$motif, plotCurve=FALSE, maxRank=5000-20, genesFormat="incidMatrix", method="aprox")$incidMatrix
  
  # collapse from motifs to tfs
  tfs = str_extract(sig_motifs$TF_highConf, "^[[:alnum:]]+")
  tfs_unique = unique(tfs)
  tfs_unique = sort(tfs_unique[!is.na(tfs_unique)])
  
  inc_mat_collapse = data.frame(matrix(ncol=dim(inc_mat)[2], nrow=0))
  colnames(inc_mat_collapse) <- colnames(inc_mat)
  
  tf_df = data.frame(tf=c(tfs_unique, "All"), frac_binding_genes=NA, mean_nes=NA, mean_auc=NA)
  
  for(i in 1:length(tfs_unique)) {
    
    t_ix = which(tfs==tfs_unique[i])
    tf_df$mean_nes[i] = mean(sig_motifs$NES[t_ix])
    tf_df$mean_auc[i] = mean(sig_motifs$AUC[t_ix])
    tf_df$frac_binding_genes[i] = length(unique(unlist(sapply(sig_motifs$enrichedGenes[t_ix], function(x) str_split(x, ";"))))) / length(gene_lists[[j]])
    
    if(length(t_ix) > 1) {
      to_add = apply(inc_mat[t_ix,], 2, function(x) round(sum(x)/length(x)))
      inc_mat_collapse = rbind(inc_mat_collapse, t(data.frame(to_add)))
    } else {
      inc_mat_collapse = rbind(inc_mat_collapse, t(as.data.frame(inc_mat[t_ix,])))
    }
  }
  
  rownames(inc_mat_collapse) = tfs_unique
  tf_df$frac_binding_genes[tf_df$tf=="All"] = sum(apply(inc_mat_collapse, 2, function(x) any(x==1))) / length(gene_lists[[j]])
  
  # construct network  
  edges = melt(as.matrix(inc_mat_collapse))
  edges = edges[which(edges[,3]==1),1:2]
  colnames(edges) = c("from","to")
  tfs = unique(as.character(edges[,1]))
  genes = unique(as.character(edges[,2]))
  nodes = data.frame(id=c(tfs, genes), label=c(tfs, genes), title=c(tfs, genes), shape=c(rep("diamond", length(tfs)), rep("elypse", length(genes))), color=c(rep("purple", length(tfs)), rep("skyblue", length(genes))))
  nodes = nodes[!duplicated(nodes$id),]
  
  # visualise
  # visNetwork(nodes, edges, main=conds[j]) %>% visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>% visPhysics(enabled=FALSE)
  
  # populate tf information for tfs_list
  tfs_list[[j]] = tf_df
  
}

all_data = lapply(tfs_list, function(x) as.character(x$tf))
# all_data = gene_lists
all_data_tbl = as.data.frame.matrix((table(stack(all_data))))
all_data_tbl = cbind(rownames(all_data_tbl), all_data_tbl)
rownames(all_data_tbl) = NULL
upset(all_data_tbl, order.by="freq", nsets=20)

# \ CONDITION ENRICHMENT --------------------------------------------------


res_cell_lines = gather(res_cell_lines, "group","score", 2:dim(res_cell_lines)[2])
res_cell_lines$type = str_extract(res_cell_lines$group, "[a-z]+$")
res_cell_lines$group = str_replace(res_cell_lines$group, "_[a-z]+$", "")

res_all = rbind(res_cell_lines, res_novakovic) # fc and p values
# split fc and p values to 2 data frames
res_all_p = filter(res_all, type=="p") %>% dplyr::select(-type) %>% spread(group, score)
res_all_fc = filter(res_all, type=="fc") %>% dplyr::select(-type) %>% spread(group, score)
rownames(res_all_p) = res_all_p$genes
rownames(res_all_fc) = res_all_fc$genes
res_all_p = res_all_p[,-1]
res_all_fc = res_all_fc[,-1]
all(names(res_all_p)==names(res_all_fc))


# CORRELATION BETWEEN CELL LINES AND PRIMARY ------------------------------

# store correlation between thp and novakovic
trts_cell_lines = c(trts_thp1, trts_u937, "5days_BG") # add 5 days bg as control
res_template = matrix(NA, nrow=length(trts_cell_lines), ncol=length(trts_novakovic))
colnames(res_template) = trts_novakovic
rownames(res_template) = trts_cell_lines

# store the results for each pathway
out_comp = vector("list", length(my_p_filter))
names(out_comp) = names(my_p_filter)

out_activity = matrix(NA, nrow=length(my_p_filter), ncol=length(c(trts_cell_lines, trts_novakovic)))
colnames(out_activity) = c(trts_cell_lines, trts_novakovic)
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
    for(k in 1:length(trts_novakovic)) {
      
      x = pathway_fc[,which(names(pathway_fc)==trts_cell_lines[j])]
      y = pathway_fc[,which(names(pathway_fc)==trts_novakovic[k])]
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
my_labels = as.character(outer(trts_cell_lines, trts_novakovic, paste, sep="_VS_"))
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

design(dds_novakovic)
dds_novakovic = DESeq(dds_novakovic, test="LRT", reduced=~donor)
res_novakovic = tbl_df(results(dds_novakovic))
res_novakovic$gene = rownames(res_novakovic)
res_novakovic$symbol = gene_list_all$hgnc_symbol
arrange(res_novakovic, padj) %>% filter(padj<1e-8)

# take the top 1000 hits
genes_top = head(res_novakovic$gene[order(res_novakovic$padj)],1e3)

for(j in 1:3) {
  y = plotCounts(dds_novakovic, gene=head(order(res_novakovic$padj),10)[j], intgroup=c("treatment","time","donor"), returnData=TRUE)
  print(ggplot(y, aes(x=time, y=count)) + geom_point(shape=17, size=3) + theme_thesis(20) + xlab("") + ylab("Count") + ggtitle(res_novakovic$gene[head(order(res_novakovic$padj),10)[j]]))
}

y = t(assays(rld_novakovic)[[1]])
dim(y)
y = y[,colnames(y) %in% genes_top]
y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
yy = data.frame(pca_res$x[,1:2])
names(yy) = c("x","y")
yy$annot_1 = paste(rld_novakovic$treatment, rld_novakovic$time, sep="_")
yy$annot_2 = rld_novakovic$donor
ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")

for_heatmap = assays(rld_novakovic)[[1]]
colnames(for_heatmap) = paste(col_data_filt$treatment, col_data_filt$time, col_data_filt$donor, sep="_")
for_heatmap = as.data.frame(for_heatmap[match(genes_top,rownames(dds_novakovic)),])
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

for_heatmap = assays(rld_novakovic)[[1]]
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

