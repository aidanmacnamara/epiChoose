
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
require(fgsea)
require(enrichR)

load("data/dat_all.RData")
load("tmp/dds_all.RData")
load("data/roi_tss.RData")


# GENE SETS ---------------------------------------------------------------

my_p = list(
  go_bp = gmtPathways("tmp/c5.bp.v6.2.symbols.gmt"), # go biological processes,
  reactome = gmtPathways("tmp/c2.cp.reactome.v6.2.symbols.gmt") # reactome
)


# GET FOLD CHANGES --------------------------------------------------------

trts = list(
  primary = c("primary_macrophage","primary_macrophage_inflamm"),
  cell_lines = c("THP-1_LPS","THP-1_PMA","THP-1_PMA+LPS","THP-1_VD3","THP-1_VD3+LPS","U937_LPS","U937_PMA","U937_PMA+LPS","U937_VD3","U937_VD3+LPS")
)

fc_res = data.frame(genes=roi_tss$hgnc_symbol)
fc_res_long = data.frame(trmt=NULL, fc=NULL, p=NULL, group=NULL)

for(i in 1:length(trts)) {
  for(j in 1:length(trts[[i]])) {
    
    if(i==1) {
      res = results(dds_all, contrast=c("group",trts[[i]][j],"primary_monocyte"))
    }
    if(i==2) {
      if(grepl("THP-1", trts[[i]][j])) {
        res = results(dds_all, contrast=c("group",trts[[i]][j],"THP-1_Baseline"))
      } else {
        res = results(dds_all, contrast=c("group",trts[[i]][j],"U937_Baseline"))
      }
    }
    
    fc_res = cbind(fc_res, res$log2FoldChange, res$padj)
    fc_res_long = rbind(fc_res_long, data.frame(trmt=trts[[i]][j], fc=res$log2FoldChange, p=res$padj, group=names(trts)[i]))
  }
}

names(fc_res)[-1] = paste(rep(unlist(trts), each=2), c("fc","p"), sep="_")
fc_res = tbl_df(fc_res)

fc_res_long$gene = rep(roi_tss$hgnc_symbol, length(unlist(trts)))
fc_res_long$trmt = as.character(fc_res_long$trmt)
fc_res_long = tbl_df(fc_res_long)

fc_res_long_filt = fc_res_long %>% filter(abs(fc) >= 1.5, p <= 0.05)
table(fc_res_long_filt$trmt)

# write_csv(fc_res_long_filt, "~/Dropbox/OTAR020/figures_dat/fc_res_long_filt.csv")
save(fc_res_long_filt, file="tmp/fc_res_long_filt.RData")

# PILE-UP PLOT ------------------------------------------------------------

# get gene list (definition of dynamic genes)
# define window
# pull in signal
# sample comps comes from fc_res_long

round_up <- function(x, to=10) {
  to*(x %/% to + as.logical(x%%to))
}

options(stringsAsFactors=FALSE)

comps = list(
  
  comp_1 = list( # everything
    sample_comp = c("THP-1_PMA","THP-1_PMA+LPS","U937_PMA","U937_PMA+LPS","primary_macrophage","primary_macrophage_inflamm"),
    fc = 1.5,
    sample_dat = data.frame(
      sample_name = c("THP-1_BR1_Baseline","THP-1_BR2_Baseline",
                      "THP-1_BR1_PMA","THP-1_BR2_PMA",
                      "THP-1_BR1_PMA+LPS","THP-1_BR2_PMA+LPS",
                      "U937_BR1_Baseline","U937_BR2_Baseline",
                      "U937_BR1_LPS","U937_BR2_LPS",
                      "U937_BR1_PMA+LPS","U937_BR2_PMA+LPS",
                      "C000S5_CD14-positive, CD16-negative classical monocyte","C0010K_CD14-positive, CD16-negative classical monocyte","C0011I_CD14-positive, CD16-negative classical monocyte","C00408_CD14-positive, CD16-negative classical monocyte","C004SQ_CD14-positive, CD16-negative classical monocyte",
                      "S0022I_macrophage","S00390_macrophage",
                      "S001MJ_inflammatory macrophage","S0022I_inflammatory macrophage"
      ),
      sample_condition = c(
        "THP-1 Baseline","THP-1 Baseline",
        "THP-1 PMA","THP-1 PMA",
        "THP-1 PMA + LPS","THP-1 PMA + LPS",
        "U937 Baseline","U937 Baseline",
        "U937 PMA","U937 PMA",
        "U937 PMA + LPS","U937 PMA + LPS",
        "Monocyte","Monocyte","Monocyte","Monocyte","Monocyte",
        "Macrophage","Macrophage",
        "Inflammatory Macrophage","Inflammatory Macrophage"
      ),
      sample_donor = NA
    )
  ),
  
  comp_2 = list( # everything
    sample_comp = c("THP-1_PMA","THP-1_PMA+LPS","U937_PMA","U937_PMA+LPS","primary_macrophage","primary_macrophage_inflamm"),
    fc = -1.5,
    sample_dat = data.frame(
      sample_name = c("THP-1_BR1_Baseline","THP-1_BR2_Baseline",
                      "THP-1_BR1_PMA","THP-1_BR2_PMA",
                      "THP-1_BR1_PMA+LPS","THP-1_BR2_PMA+LPS",
                      "U937_BR1_Baseline","U937_BR2_Baseline",
                      "U937_BR1_LPS","U937_BR2_LPS",
                      "U937_BR1_PMA+LPS","U937_BR2_PMA+LPS",
                      "C000S5_CD14-positive, CD16-negative classical monocyte","C0010K_CD14-positive, CD16-negative classical monocyte","C0011I_CD14-positive, CD16-negative classical monocyte","C00408_CD14-positive, CD16-negative classical monocyte","C004SQ_CD14-positive, CD16-negative classical monocyte",
                      "S0022I_macrophage","S00390_macrophage",
                      "S001MJ_inflammatory macrophage","S0022I_inflammatory macrophage"
      ),
      sample_condition = c(
        "THP-1 Baseline","THP-1 Baseline",
        "THP-1 PMA","THP-1 PMA",
        "THP-1 PMA + LPS","THP-1 PMA + LPS",
        "U937 Baseline","U937 Baseline",
        "U937 PMA","U937 PMA",
        "U937 PMA + LPS","U937 PMA + LPS",
        "Monocyte","Monocyte","Monocyte","Monocyte","Monocyte",
        "Macrophage","Macrophage",
        "Inflammatory Macrophage","Inflammatory Macrophage"
      ),
      sample_donor = NA
    )
  )
  
)

# store the results of all the comparisons in these objects
total_signal = data.frame()
total_diff = data.frame()
total_gene_orders = vector("list", length(comps))
names(total_gene_orders) = names(comps)

for(i in 1:length(comps)) {
  
  print(paste("Comparison", i))
  
  # get the bigwig files
  files = dat_all$tss$H3K27ac$annot$Bigwig[match(comps[[i]]$sample_dat$sample_name, dat_all$tss$H3K27ac$annot$Label)]
  files = str_replace(files, "/GWD/bioinfo/projects/", "/Volumes/am673712/links/")
  
  # define dynamic genes for comparisons
  dynamic_genes = vector("list",length(comps[[i]]$sample_comp))
  names(dynamic_genes) = comps[[i]]$sample_comp
  
  # get any genes with fc > 1.5 across all the comparisons for each plot
  for(j in 1:length(comps[[i]]$sample_comp)) {
    if(comps[[i]]$fc > 0) {
      dynamic_genes[[j]] = fc_res_long %>% filter(p <= 0.05, fc >= comps[[i]]$fc, trmt==comps[[i]]$sample_comp[j]) %>% dplyr::select(gene) %>% unlist() %>% as.character()
    } else {
      dynamic_genes[[j]] = fc_res_long %>% filter(p <= 0.05, fc <= comps[[i]]$fc, trmt==comps[[i]]$sample_comp[j]) %>% dplyr::select(gene) %>% unlist() %>% as.character()
    }
  }
  
  # plot the venn diagram
  plot(euler(dynamic_genes), quantities=TRUE)
  
  # get these genes into long format 'gene_tbl_long'
  gene_tbl = as.data.frame.matrix((table(stack(dynamic_genes))))
  gene_tbl$gene = rownames(gene_tbl); gene_tbl = tbl_df(gene_tbl)
  gene_tbl_long = gene_tbl %>% gather("sample","deregulated",1:length(comps[[i]]$sample_comp))
  gene_tbl_long$comp = names(comps)[i]
  
  # match genes to granges total genome (roi_tss)
  dynamic_genes = unique(unlist(dynamic_genes))
  gene_list = roi_tss[match(dynamic_genes,roi_tss$hgnc_symbol)]
  
  # increase the window
  start(gene_list) = start(gene_list) - 1e4
  end(gene_list) = end(gene_list) + 1e4
  
  # also store the list of genes for plot levels
  total_gene_orders[[i]] = dynamic_genes
  
  # define the number of tiles per gene for the pile-up plot
  n_tiles = 21
  
  # grab the signal values
  y = bplapply(1:length(files), function(x) rtracklayer::import(con=files[x], which=gene_list))
  names(y) = comps[[i]]$sample_dat$sample_name
  
  # tile each gene
  gene_list_tiled = tile(gene_list, n=n_tiles)
  
  for(j in 1:length(gene_list_tiled)) {
    gene_list_tiled[[j]]$gene = gene_list$hgnc_symbol[j]
  }
  gene_list_tiled = unlist(gene_list_tiled)
  mcol_ix = 2
  
  for(j in 1:length(y)) { # for each gene, split the signal across the tiles
    
    print(paste("sample:", j))
    
    ols = findOverlaps(gene_list_tiled, y[[j]])
    binned_score = rep(NA, length(gene_list_tiled))
    for(k in 1:length(binned_score)) {
      binned_score[k] = mean(y[[j]]$score[subjectHits(ols)[queryHits(ols)==k]], na.rm=TRUE)
    }
    binned_score[is.nan(binned_score)] = 0
    mcols(gene_list_tiled)[,mcol_ix] = binned_score
    names(mcols(gene_list_tiled))[mcol_ix] = comps[[i]]$sample_dat$sample_name[j]
    mcol_ix = mcol_ix+1
  }
  
  tile_size = round_up((width(gene_list[1])/n_tiles), to=100)
  coord_labels = c(paste0("-", rev(tile_size*1:((n_tiles-1)/2))),"TSS",paste0("+", tile_size*1:((n_tiles-1)/2)))
  gene_list_tiled$coord_ix = rep(1:n_tiles, length(gene_list))
  
  # make the df from the granges object
  
  gene_list_df = tbl_df(as.data.frame(gene_list_tiled))
  names(gene_list_df)[7:(dim(gene_list_df)[2]-1)] = as.character(comps[[i]]$sample_dat$sample_name) # correct name
  
  # change to long format
  
  gene_list_df_long = gene_list_df %>% gather("sample","score",7:(dim(gene_list_df)[2]-1))
  gene_list_df_long$comp = names(comps)[i]
  
  # collapse replicates
  
  gene_list_df_long$sample_condition = comps[[i]]$sample_dat$sample_condition[match(gene_list_df_long$sample, comps[[i]]$sample_dat$sample_name)]
  gene_list_df_long_summ = gene_list_df_long %>% group_by(seqnames,start,end,width,strand,gene,coord_ix,sample_condition,comp) %>% summarise(score_mean=mean(score))
  
  # add to totals
  total_signal = rbind(total_signal, as.data.frame(gene_list_df_long_summ))
  total_diff = rbind(total_diff, gene_tbl_long)
  
}

save(total_signal, file="tmp/total_signal.RData") # savepoint
save(total_diff, file="tmp/total_diff.RData") # savepoint
save(total_gene_orders, file="tmp/total_gene_orders.RData") # savepoint
save(comps, file="tmp/comps.RData")


# GENE ENRICHMENT ---------------------------------------------------------

load("tmp/total_diff.RData")
load("tmp/total_signal.RData")
load("tmp/total_gene_orders.RData")
load("tmp/fc_res_long_filt.RData")

table(fc_res_long_filt$trmt)
fc_res_long_filt$dir = ifelse(fc_res_long_filt$fc > 0, "up", "down")
fc_res_long_filt$trmt_dir = paste(fc_res_long_filt$trmt, fc_res_long_filt$dir, sep="_")

gene_sets = sapply(unique(fc_res_long_filt$trmt_dir), function(x) fc_res_long_filt$gene[fc_res_long_filt$trmt_dir==x])
lapply(gene_sets, length)

# enrichment

dbs = tbl_df(listEnrichrDbs())
dbs = c("GO_Biological_Process_2017","Reactome_2016")

enrich_process <- function(x) {
  res = enrichr(unique(x), dbs)
  res = lapply(res, function(x) tbl_df(filter(x, `Adjusted.P.value` <= 0.05)))
}

enrich_res <- lapply(gene_sets, enrich_process)

# gene set enrichment analysis

gsea_res = vector("list", length(gene_sets))
names(gsea_res) = names(gene_sets)

for(i in 1:length(gene_sets)) {
  
  gsea_dat = fc_res_long_filt %>% filter(trmt_dir==names(gene_sets)[i]) %>% select(gene, fc) %>% mutate(fc=abs(fc)) %>% arrange(desc(fc))
  
  gsea_in = gsea_dat$fc
  names(gsea_in) = gsea_dat$gene
  
  if(length(gsea_in) < 2) next
  
  res_gsea <- lapply(my_p, function(x) fgsea(x, stats=gsea_in, nperm=1000))
  gsea_res[[i]] = lapply(res_gsea, function(x) x %>% as_tibble() %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj))
  
}


# BIMODAL SELECTION -------------------------------------------------------

require(mclust)
require(cowplot) 

# poc

my_sample = "THP-1_BR1_Baseline"
dat = dat_all$tss$H3K27ac$res

to_plot = data.frame(
  gene = roi_tss$hgnc_symbol,
  signal = log(dat[which(rownames(dat)==my_sample),]),
  expression = log(unlist(dat_all$tss$RNA$res[which(rownames(dat_all$tss$RNA$res)==my_sample),]))
)

to_plot = to_plot[!is.infinite(to_plot$signal),]

gmm = Mclust(to_plot$signal, G=2)
to_plot$class = factor(gmm$classification, labels=c("High","Low"))

# main plot if primary (no rna)

ggplot() + geom_density(data=to_plot, aes(x=signal, fill=class), alpha=0.2, size=0.2) + ggpubr::fill_palette("jco")

# main plot if cell line

p_main = ggplot(to_plot, aes(x=signal, y=expression, color=class)) + geom_point(alpha=0.2) + theme_thesis(20) + ggpubr::color_palette("jco")
x_dens = axis_canvas(p_main, axis="x") + geom_density(data=to_plot, aes(x=signal, fill=class), alpha=0.2, size=0.2) + ggpubr::fill_palette("jco") # marginal densities along x-axis
y_dens <- axis_canvas(p_main, axis="y", coord_flip=TRUE) + geom_density(data=to_plot, aes(x=expression, fill=class), alpha=0.2, size=0.2) + coord_flip() + ggpubr::fill_palette("jco") # marginal densities along y-axis
p_1 = insert_xaxis_grob(p_main, x_dens, grid::unit(.2, "null"), position="top")
p_2 = insert_yaxis_grob(p_1, y_dens, grid::unit(.2, "null"), position="right")
ggdraw(p_2)

# across the different groups - define on/off genes

my_groups = unique(as.character(comps$comp_1$sample_dat$sample_condition))
my_groups_data = vector("list", length(my_groups))
names(my_groups_data) = my_groups

dat = dat_all$tss$H3K27ac$res

for(i in 1:length(my_groups)) {
  
  samples = comps$comp_1$sample_dat$sample_name[comps$comp_1$sample_dat$sample_condition==my_groups[i]]
  to_plot = data.frame()
  
  for(j in 1:length(samples)) {
    
    to_plot_sample = data.frame(
      gene = roi_tss$hgnc_symbol,
      signal = log(dat[which(rownames(dat)==samples[j]),]),
      expression = log(unlist(dat_all$tss$RNA$res[which(rownames(dat_all$tss$RNA$res)==samples[j]),])),
      sample = samples[j]
    )
    
    to_plot_sample = tbl_df(to_plot_sample[!is.infinite(to_plot_sample$signal),])
    
    gmm = Mclust(to_plot_sample$signal, G=2)
    to_plot_sample$class = factor(gmm$classification, labels=c("High","Low"))
    
    to_plot = rbind(to_plot, to_plot_sample)
  }
  
  high_group = to_plot %>% group_by(gene) %>% summarise(high = all(class=="High"))
  high_group = high_group %>% filter(high==TRUE) %>% select(gene) %>% unlist() %>% as.character()
  high_group = high_group[high_group!=""]
  
  low_group = to_plot %>% group_by(gene) %>% summarise(low = all(class=="Low"))
  low_group = low_group %>% filter(low==TRUE) %>% select(gene) %>% unlist() %>% as.character()
  low_group = low_group[low_group!=""]
  
  my_groups_data[[i]] = list(total=to_plot, high_group=high_group, low_group=low_group)
  
}

# how many high signal per sample
lapply(my_groups_data, function(x) x[[1]] %>% filter(class=="High") %>% group_by(sample) %>% summarise(N=n()))

# how many high signal per group
lapply(my_groups_data, function(x) length(x$high_group))

# how many low signal per group
lapply(my_groups_data, function(x) length(x$low_group))

# plot above
plot(euler(lapply(my_groups_data, function(x) x$high_group)))

# check enrichment of primary groups

enrich_process <- function(x) {
  res = enrichr(unique(x), dbs)
  res = lapply(res, function(x) tbl_df(filter(x, `Adjusted.P.value` <= 0.05)))
  # res = lapply(res, function(x) tbl_df(head(x)))
}

enrich_res <- lapply(lapply(my_groups_data[7:9], function(x) x$high_group), enrich_process)

# send list to team

gene_intersects = data.frame()

primary_high = Reduce(intersect, lapply(my_groups_data[7:9], function(x) x$high_group))
primary_low = Reduce(intersect, lapply(my_groups_data[7:9], function(x) x$low_group))
thp_high = Reduce(intersect, lapply(my_groups_data[1:3], function(x) x$high_group))
thp_low = Reduce(intersect, lapply(my_groups_data[1:3], function(x) x$low_group))

overlaps_1 = list(
  pma_up = fc_res_long_filt %>% filter(trmt=="THP-1_PMA", fc>0) %>% select(gene) %>% unlist %>% as.character(),
  primary_high = primary_high
)
plot(euler(overlaps_1))
gene_intersects = rbind(gene_intersects, data.frame(comp="pma_up_vs_primary_high", genes=Reduce(intersect, overlaps_1)))

overlaps_2 = list(
  pma_down = fc_res_long_filt %>% filter(trmt=="THP-1_PMA", fc<0) %>% select(gene) %>% unlist %>% as.character(),
  primary_low = primary_low
)
plot(euler(overlaps_2))
gene_intersects = rbind(gene_intersects, data.frame(comp="pma_down_vs_primary_low", genes=Reduce(intersect, overlaps_2)))

overlaps_3 = list(
  primary_up = fc_res_long_filt %>% filter(trmt=="primary_macrophage", fc>0) %>% select(gene) %>% unlist %>% as.character(),
  thp_high = thp_high
)
plot(euler(overlaps_3))
gene_intersects = rbind(gene_intersects, data.frame(comp="primary_up_vs_thp_high", genes=Reduce(intersect, overlaps_3)))

overlaps_4 = list(
  primary_down = fc_res_long_filt %>% filter(trmt=="primary_macrophage", fc<0) %>% select(gene) %>% unlist %>% as.character(),
  thp_low = thp_low
)
plot(euler(overlaps_4))
gene_intersects = rbind(gene_intersects, data.frame(comp="primary_down_vs_thp_low", genes=Reduce(intersect, overlaps_4)))

table(gene_intersects$comp)
write_csv(gene_intersects, "~/Dropbox/OTAR020/figures_dat/gene_intersects.csv")


# CHECK AGAINST PAPER -----------------------------------------------------

promoter_peaks = read_excel("tmp/1-s2.0-S0092867416313162-mmc1.xlsx", sheet="S1B. H3K27ac AcP", skip=3)


# TF ANALYSIS -------------------------------------------------------------

tfs_list =  vector("list", length(gene_sets))
names(tfs_list) = names(gene_sets)

# import motif-gene rankings
# how is this ranking performed?
# how much do the results change varying the tss window?
motif_rankings <- importRankings("tmp/hg19-tss-centered-5kb-7species.mc9nr.feather")

# import motif-tf annotations
data(motifAnnotations_hgnc)

# run enrichment
# question: what motifs are enriched in each cluster?
motif_enrichment <- cisTarget(gene_sets, motif_rankings, motifAnnot=motifAnnotations_hgnc)

for(j in 1:length(gene_sets)) { # visualise
  
  print(j)
  
  # get the significant motifs
  # can also apply motif-tf  mapping threshold here
  sig_motifs = filter(motif_enrichment, geneSet==names(gene_sets)[j], NES>=4)
  
  if(dim(sig_motifs)[1]==0) next
  
  # get all interactions between significant tfs and regulated genes
  # this is also identifying the significant regulated genes
  inc_mat = getSignificantGenes(gene_sets[[j]], motif_rankings, signifRankingNames=sig_motifs$motif, plotCurve=FALSE, maxRank=5000-20, genesFormat="incidMatrix", method="aprox")$incidMatrix
  
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
    tf_df$frac_binding_genes[i] = length(unique(unlist(sapply(sig_motifs$enrichedGenes[t_ix], function(x) str_split(x, ";"))))) / length(gene_sets[[j]])
    
    if(length(t_ix) > 1) {
      to_add = apply(inc_mat[t_ix,], 2, function(x) round(sum(x)/length(x)))
      inc_mat_collapse = rbind(inc_mat_collapse, t(data.frame(to_add)))
    } else {
      inc_mat_collapse = rbind(inc_mat_collapse, t(as.data.frame(inc_mat[t_ix,])))
    }
  }
  
  rownames(inc_mat_collapse) = tfs_unique
  tf_df$frac_binding_genes[tf_df$tf=="All"] = sum(apply(inc_mat_collapse, 2, function(x) any(x==1))) / length(gene_sets[[j]])
  
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


# THP-1 VS PMA CHARACTERISATION -------------------------------------------

my_p = gmtPathways("tmp/c2.cp.reactome.v6.2.symbols.gmt") # biological processes
# my_p = gmtPathways("tmp/c5.bp.v6.2.symbols.gmt")

my_contrasts = list(
  c("THP-1_PMA","THP-1_Baseline"),
  c("THP-1_VD3","THP-1_Baseline"),
  c("U937_PMA","U937_Baseline"),
  c("U937_VD3","U937_Baseline"),
  c("THP-1_PMA","U937_PMA")
)

sig_pathways = vector("list", length(my_contrasts)*2)
names(sig_pathways) = paste(c(rep("chip",length(my_contrasts)),rep("rna",length(my_contrasts))), unlist(lapply(my_contrasts, function(x) paste(x, collapse="_vs_"))), sep="_")
sig_diff = sig_pathways
sig_count = 1

for(j in 1:2) {
  
  if(j==1) dds_to_use = dds_cell_lines
  if(j==2) dds_to_use = dds_cell_lines_rna
  
  for(i in 1:length(my_contrasts)) {
    
    res = results(dds_to_use, contrast=c("cell_condition", my_contrasts[[i]]))
    
    if(j==1) res$gene = gene_list_all$hgnc_symbol
    if(j==2) res$gene = gene_list_all$hgnc_symbol[-gene_missing_ix]
    
    res = tbl_df(res) %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
    
    gsea_in = res$log2FoldChange
    names(gsea_in) = res$gene
    res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
    res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES)) %>% filter(padj < 0.1)
    res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()
    
    sig_pathways[[sig_count]] = res_tidy$pathway
    if(j==1) sig_diff[[sig_count]] = filter(res, log2FoldChange>1) %>% select(gene) %>% unlist() %>% as.character()
    if(j==2) sig_diff[[sig_count]] = filter(res, log2FoldChange>3) %>% select(gene) %>% unlist() %>% as.character()
    sig_count = sig_count + 1
  }
}

lapply(sig_pathways, length)
lapply(sig_diff, length)


# PLOT --------------------------------------------------------------------

plot(euler(sig_pathways), quantities=TRUE)
plot(euler(sig_diff), quantities=TRUE)

my_dat_genes_all_tbl = as.data.frame.matrix((table(stack(sig_diff))))
my_dat_genes_all_tbl = cbind(rownames(my_dat_genes_all_tbl), my_dat_genes_all_tbl)
rownames(my_dat_genes_all_tbl) = NULL
upset(my_dat_genes_all_tbl, order.by="freq", nsets=10)

x_ix = which(colData(dds_cell_lines)$Label=="THP-1_BR1_Baseline")
y_ix = which(colData(dds_cell_lines)$Label=="U937_BR1_VD3")

data.frame(chip=assays(rld_cell_lines)[[1]][-gene_missing_ix,x_ix], rna=assays(rld_cell_lines_rna)[[1]][,y_ix]) %>% ggplot(aes(x=chip,y=rna)) + geom_point(size=0.5, alpha=0.2) + theme_thesis()
