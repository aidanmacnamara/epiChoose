
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
load("tmp/dds_primary.RData")

load("tmp/rld_cell_lines.RData")
load("tmp/rld_saeed.RData")
load("tmp/rld_novakovic.RData")
load("tmp/rld_primary.RData")

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


# GET FOLD CHANGES --------------------------------------------------------

novakovic = as.character(unique(dds_novakovic$time_treatment))
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

fc_res = tbl_df(fc_res)

fc_res_long$gene = rep(gene_list_all$hgnc_symbol, length(unlist(trts)))
names(fc_res)[-1] = paste(rep(unlist(trts), each=2), c("fc","p"), sep="_")
fc_res_long$trmt = as.character(fc_res_long$trmt)
fc_res_long = tbl_df(fc_res_long)


# DYNAMIC REGIONS ---------------------------------------------------------

# annotate with significant changes from monocytes to macrophages (late vs. early)
res_early_late = results(dds_primary, contrast=c("group","late","early"))
res_early_late$gene = gene_list_all$hgnc_symbol
res_early_late = tbl_df(res_early_late) %>% filter(padj < 0.05, abs(log2FoldChange) > 1) %>% arrange(desc(abs(log2FoldChange)))

up_primary = filter(res_early_late, log2FoldChange>0) %>% select(gene) %>% unlist %>% as.character()
down_primary = filter(res_early_late, log2FoldChange<0) %>% select(gene) %>% unlist %>% as.character()
up_cell_line_t = filter(fc_res_long, p <= 0.05, fc >= 1, trmt=="THP-1_PMA") %>% dplyr::select(gene) %>% unlist() %>% as.character()
down_cell_line_t = filter(fc_res_long, p <= 0.05, fc <= -1, trmt=="THP-1_PMA") %>% dplyr::select(gene) %>% unlist() %>% as.character()
up_cell_line_u = filter(fc_res_long, p <= 0.05, fc >= 1, trmt=="U937_PMA") %>% dplyr::select(gene) %>% unlist() %>% as.character()
down_cell_line_u = filter(fc_res_long, p <= 0.05, fc <= -1, trmt=="U937_PMA") %>% dplyr::select(gene) %>% unlist() %>% as.character()

dyn_genes = unique(
  c(up_primary,up_cell_line_t,up_cell_line_u,down_primary,down_cell_line_t,down_cell_line_u)
)


# HEATMAP -----------------------------------------------------------------

# get the log-normalised signals
y = assays(rld_all)[[1]]
dim(y)

# define indices
gene_ix = match(dyn_genes, gene_list_all$hgnc_symbol)
rownames(y) = gene_list_all$hgnc_symbol
y = y[gene_ix,]
dim(y)
y = y[,!grepl("BG|LPS|glucan|broad", colnames(y), ignore.case=TRUE)]
dim(y)

# make annotation for heatmap
annot_bar = data.frame(name=rownames(y), primary_down=0, primary_up=0, thp_down=0, thp_up=0, u937_down=0, u937_up=0)
annot_bar$primary_down[annot_bar$name %in% down_primary] = 1
annot_bar$thp_down[annot_bar$name %in% down_cell_line_t] = 1
annot_bar$u937_down[annot_bar$name %in% down_cell_line_u] = 1
annot_bar$primary_up[annot_bar$name %in% up_primary] = 1
annot_bar$thp_up[annot_bar$name %in% up_cell_line_t] = 1
annot_bar$u937_up[annot_bar$name %in% up_cell_line_u] = 1

pheatmap(t(y), annotation_col=annot_bar, show_colnames=FALSE)

# plot
heatmaply(t(y), col=bluered(75), cexCol=0.5, col_side_colors=annot_bar)

# correlation matrix
yy = cor(y)
heatmaply(yy, col=bluered(75))

# look at primary data only
row_ix = grep("SANQUIN|N000", rownames(yy))
yyy = yy[row_ix, -row_ix]


# TF ANALYSIS -------------------------------------------------------------

fc_res_out = filter(fc_res, genes %in% dyn_genes) %>% dplyr::select(-grep("BG|LPS|VD3|glucan", names(.)))
fc_res_out = cbind(fc_res_out, primary_down=0, primary_up=0, thp_down=0, thp_up=0, u937_down=0, u937_up=0)
fc_res_out$primary_down[fc_res_out$genes %in% down_primary] = 1
fc_res_out$thp_down[fc_res_out$genes %in% down_cell_line_t] = 1
fc_res_out$u937_down[fc_res_out$genes %in% down_cell_line_u] = 1
fc_res_out$primary_up[fc_res_out$genes %in% up_primary] = 1
fc_res_out$thp_up[fc_res_out$genes %in% up_cell_line_t] = 1
fc_res_out$u937_up[fc_res_out$genes %in% up_cell_line_u] = 1

plot(euler(fc_res_out[,14:19]), quantities=TRUE)

# what is the enrichment for a condition, e.g. thp1 + pma
conds = c("Primary", "THP-1", "U937")
conds_list = vector("list", length(conds))
names(conds_list) = conds

tfs_list =  vector("list", length(conds))
names(tfs_list) = conds

conds_list$Primary = res_early_late %>% filter(padj <= 0.05, log2FoldChange >= 2) %>% dplyr::select(gene) %>% unlist() %>% as.character()
conds_list$`THP-1` = fc_res_long %>% filter(trmt=="THP-1_PMA", fc >= 2, p <= 0.05) %>% dplyr::select(gene) %>% unlist() %>% as.character()
conds_list$U937 = fc_res_long %>% filter(trmt=="U937_PMA", fc >= 2, p <= 0.05) %>% dplyr::select(gene) %>% unlist() %>% as.character()

# import motif-gene rankings
# how is this ranking performed?
# how much do the results change varying the tss window?
motif_rankings <- importRankings("tmp/hg19-tss-centered-5kb-7species.mc9nr.feather")

# import motif-tf annotations
data(motifAnnotations_hgnc)

# run enrichment
# question: what motifs are enriched in each cluster?
motif_enrichment <- cisTarget(conds_list, motif_rankings, motifAnnot=motifAnnotations_hgnc)

for(j in 1:length(conds_list)) { # visualise
  
  print(j)
  
  # get the significant motifs
  # can also apply motif-tf  mapping threshold here
  sig_motifs = filter(motif_enrichment, geneSet==names(conds_list)[j], NES>=4)
  
  if(dim(sig_motifs)[1]==0) next
  
  # get all interactions between significant tfs and regulated genes
  # this is also identifying the significant regulated genes
  inc_mat = getSignificantGenes(conds_list[[j]], motif_rankings, signifRankingNames=sig_motifs$motif, plotCurve=FALSE, maxRank=5000-20, genesFormat="incidMatrix", method="aprox")$incidMatrix
  
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
    tf_df$frac_binding_genes[i] = length(unique(unlist(sapply(sig_motifs$enrichedGenes[t_ix], function(x) str_split(x, ";"))))) / length(conds_list[[j]])
    
    if(length(t_ix) > 1) {
      to_add = apply(inc_mat[t_ix,], 2, function(x) round(sum(x)/length(x)))
      inc_mat_collapse = rbind(inc_mat_collapse, t(data.frame(to_add)))
    } else {
      inc_mat_collapse = rbind(inc_mat_collapse, t(as.data.frame(inc_mat[t_ix,])))
    }
  }
  
  rownames(inc_mat_collapse) = tfs_unique
  tf_df$frac_binding_genes[tf_df$tf=="All"] = sum(apply(inc_mat_collapse, 2, function(x) any(x==1))) / length(conds_list[[j]])
  
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


# PATHWAY COMPARISON ------------------------------------------------------

# split fc and p values to 2 data frames
fc_res_p = fc_res %>% dplyr::select(matches("_p$")) %>% as.data.frame()
names(fc_res_p) = str_replace(names(fc_res_p), "_p$", "")
rownames(fc_res_p) = fc_res$genes

fc_res_fc = fc_res %>% dplyr::select(matches("_fc$")) %>% as.data.frame()
names(fc_res_fc) = str_replace(names(fc_res_fc), "_fc$", "")
rownames(fc_res_fc) = fc_res$genes

all(names(fc_res_p)==names(fc_res_fc))

# store correlation between thp and novakovic
trts_cell_lines = c(trts$cell_lines, "5days_BG") # add 5 days bg as control
res_template = matrix(NA, nrow=length(trts_cell_lines), ncol=length(trts$novakovic))
colnames(res_template) = trts$novakovic
rownames(res_template) = trts_cell_lines

# store the results for each pathway
out_comp = vector("list", length(my_p_filter))
names(out_comp) = names(my_p_filter)

out_activity = matrix(NA, nrow=length(my_p_filter), ncol=length(c(trts_cell_lines, trts$novakovic)))
colnames(out_activity) = c(trts_cell_lines, trts$novakovic)
rownames(out_activity) = names(my_p_filter)

for(i in 1:length(out_comp)) {
  
  res_copy = res_template
  g_ix = which(rownames(fc_res_fc) %in% my_p_filter[[i]]) # get the gene indices for the pathway
  if(is_empty(g_ix)) next
  
  pathway_fc = fc_res_fc[g_ix,] # get the fold changes
  pathway_p = fc_res_p[g_ix,] # get the p-values
  my_p_annot$mean_fc[i] = mean(as.numeric(unlist(pathway_fc)), na.rm=TRUE) # what is the mean fc for the pathway?
  my_p_annot$sig_fcs[i] = sum(pathway_p<0.05, na.rm=TRUE) / length(unlist(pathway_p)) # how many sig fcs for the pathway?
  
  # pathway_fc[pathway_p > 0.05] = NA # set non-significant (padj) fcs to na
  
  fcs_per_condition_func <- function(x) {
    if(is.na(x)) {
      return
    }
  }
  
  fcs_per_condition = apply(pathway_p, 2, function(x) sum(x < 0.05, na.rm=TRUE) / dim(pathway_p)[1])
  out_activity[i,] = fcs_per_condition[match(colnames(out_activity),names(fcs_per_condition))]
  
  for(j in 1:length(trts_cell_lines)) {
    for(k in 1:length(trts$novakovic)) {
      
      x = pathway_fc[,which(names(pathway_fc)==trts_cell_lines[j])]
      y = pathway_fc[,which(names(pathway_fc)==trts$novakovic[k])]
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
my_labels = as.character(outer(trts_cell_lines, trts$novakovic, paste, sep="_VS_"))
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
  g_ix = which(rownames(fc_res_fc) %in% my_p_filter[[i]])
  dat_fc = fc_res_fc[g_ix,]
  pheatmap(dat_fc, main=names(my_p_filter)[i], fontsize_row=5)
  
}

gsea_in = fc_res_fc$`5days_LPS`
names(gsea_in) = rownames(fc_res_fc)
res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()


