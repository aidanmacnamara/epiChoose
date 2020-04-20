
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
require(EnhancedVolcano)

load("data/dat_all.RData")
load("tmp/dds_all.RData")
load("data/roi_tss.RData")


# GENE SETS ---------------------------------------------------------------

my_p = list(
  go_bp = gmtPathways("tmp/gene_sets/c5.bp.v6.2.symbols.gmt"), # go biological processes,
  reactome = gmtPathways("tmp/gene_sets/c2.cp.reactome.v6.2.symbols.gmt") # reactome
)


# GET FOLD CHANGES --------------------------------------------------------

trts = list(
  primary = c("primary_macrophage","primary_macrophage_inflamm"),
  cell_lines = c("THP-1_LPS","THP-1_PMA","THP-1_PMA+LPS","THP-1_VD3","THP-1_VD3+LPS","U937_LPS","U937_PMA","U937_PMA+LPS","U937_VD3","U937_VD3+LPS")
)
trts$cell_lines = str_replace_all(trts$cell_lines, "[-+]" ,"_")

fc_res = data.frame(genes=roi_tss$hgnc_symbol)
fc_res_long = data.frame(trmt=NULL, fc=NULL, p=NULL, group=NULL)

for(i in 1:length(trts)) {
  for(j in 1:length(trts[[i]])) {
    
    if(i==1) {
      contr_list = paste0("group", c(trts[[i]][j],"primary_monocyte"))
    }
    if(i==2) {
      if(grepl("THP_1", trts[[i]][j])) {
        contr_list = paste0("group", c(trts[[i]][j],"THP_1_Baseline"))
      } else {
        contr_list = paste0("group", c(trts[[i]][j],"U937_Baseline"))
      }
    }
    
    my_args = list(
      paste(contr_list[1], contr_list[2], sep=" - "),
      levels=colnames(coefficients(dds_all))
    )
    contr = do.call(makeContrasts, my_args)
    contr_fit = contrasts.fit(dds_all, contr)
    contr_fit = eBayes(contr_fit)
    
    fc_res = cbind(fc_res, contr_fit$coefficients[,1], contr_fit$F.p.value)
    fc_res_long = rbind(fc_res_long, data.frame(trmt=trts[[i]][j], fc=contr_fit$coefficients[,1], p=contr_fit$F.p.value, group=names(trts)[i]))
  }
}

names(fc_res)[-1] = paste(rep(unlist(trts), each=2), c("fc","p"), sep="_")
fc_res = tbl_df(fc_res)

fc_res_long$gene = rep(roi_tss$hgnc_symbol, length(unlist(trts)))
fc_res_long$trmt = as.character(fc_res_long$trmt)
fc_res_long = tbl_df(fc_res_long)
save(fc_res_long, file="tmp/fc_res_long.RData")

volc_plot = data.frame(filter(fc_res_long,  trmt=="primary_macrophage"))
rownames(volc_plot) = volc_plot$gene
EnhancedVolcano(volc_plot, lab=rownames(volc_plot), x="fc", y="p", FCcutoff=3, legendLabels=c("","","","FC > 3"))

fc_res_long_filt = fc_res_long %>% filter(abs(fc) >= 2, p <= 0.05)
table(fc_res_long_filt$trmt)

write_csv(fc_res_long_filt, "~/Dropbox/OTAR020/figures_dat/fc_res_long_filt.csv")
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
    sample_comp = c("THP_1_PMA","THP_1_PMA_LPS","U937_PMA","U937_PMA_LPS","primary_macrophage","primary_macrophage_inflamm"),
    fc = 1.5,
    sample_dat = data.frame(
      sample_name = c("THP-1_BR1_Baseline","THP-1_BR2_Baseline",
                      "THP-1_BR1_PMA","THP-1_BR2_PMA",
                      "THP-1_BR1_PMA+LPS","THP-1_BR2_PMA+LPS",
                      "U937_BR1_Baseline","U937_BR2_Baseline",
                      "U937_BR1_LPS","U937_BR2_LPS",
                      "U937_BR1_PMA+LPS","U937_BR2_PMA+LPS",
                      "C000S5_CD14-positive, CD16-negative classical monocyte","C0010K_CD14-positive, CD16-negative classical monocyte","C0011I_CD14-positive, CD16-negative classical monocyte","C00408_CD14-positive, CD16-negative classical monocyte","C004SQ_CD14-positive, CD16-negative classical monocyte","Monocyte_11","Monocyte_12",
                      "S0022I_macrophage","S00390_macrophage","Macrophage_11","Macrophage_16",
                      "S001MJ_inflammatory macrophage","S0022I_inflammatory macrophage"
      ),
      sample_condition = c(
        "THP-1 Baseline","THP-1 Baseline",
        "THP-1 PMA","THP-1 PMA",
        "THP-1 PMA + LPS","THP-1 PMA + LPS",
        "U937 Baseline","U937 Baseline",
        "U937 PMA","U937 PMA",
        "U937 PMA + LPS","U937 PMA + LPS",
        "Monocyte","Monocyte","Monocyte","Monocyte","Monocyte","Monocyte","Monocyte",
        "Macrophage","Macrophage","Macrophage","Macrophage",
        "Inflammatory Macrophage","Inflammatory Macrophage"
      ),
      sample_donor = NA
    )
  ),
  
  comp_2 = list( # everything
    sample_comp = c("THP_1_PMA","THP_1_PMA_LPS","U937_PMA","U937_PMA_LPS","primary_macrophage","primary_macrophage_inflamm"),
    fc = -1.5,
    sample_dat = data.frame(
      sample_name = c("THP-1_BR1_Baseline","THP-1_BR2_Baseline",
                      "THP-1_BR1_PMA","THP-1_BR2_PMA",
                      "THP-1_BR1_PMA+LPS","THP-1_BR2_PMA+LPS",
                      "U937_BR1_Baseline","U937_BR2_Baseline",
                      "U937_BR1_LPS","U937_BR2_LPS",
                      "U937_BR1_PMA+LPS","U937_BR2_PMA+LPS",
                      "C000S5_CD14-positive, CD16-negative classical monocyte","C0010K_CD14-positive, CD16-negative classical monocyte","C0011I_CD14-positive, CD16-negative classical monocyte","C00408_CD14-positive, CD16-negative classical monocyte","C004SQ_CD14-positive, CD16-negative classical monocyte","Monocyte_11","Monocyte_12",
                      "S0022I_macrophage","S00390_macrophage","Macrophage_11","Macrophage_16",
                      "S001MJ_inflammatory macrophage","S0022I_inflammatory macrophage"
      ),
      sample_condition = c(
        "THP-1 Baseline","THP-1 Baseline",
        "THP-1 PMA","THP-1 PMA",
        "THP-1 PMA + LPS","THP-1 PMA + LPS",
        "U937 Baseline","U937 Baseline",
        "U937 PMA","U937 PMA",
        "U937 PMA + LPS","U937 PMA + LPS",
        "Monocyte","Monocyte","Monocyte","Monocyte","Monocyte","Monocyte","Monocyte",
        "Macrophage","Macrophage","Macrophage","Macrophage",
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
    
    # binned_score = rep(NA, length(gene_list_tiled))
    # for(k in 1:length(binned_score)) {
    #   binned_score[k] = mean(y[[j]]$score[subjectHits(ols)[queryHits(ols)==k]], na.rm=TRUE)
    # }
    
    binned_score = bplapply(1:length(gene_list_tiled), function(k) mean(y[[j]]$score[subjectHits(ols)[queryHits(ols)==k]], na.rm=TRUE)) # parallelise to speed up
    binned_score = unlist(binned_score)
    
    binned_score[is.na(binned_score)] = 0
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
  
  q_norm = FALSE
  if(q_norm) {
    require(limma)
    dat_to_m = as.matrix(gene_list_df[,7:(dim(gene_list_df)[2]-1)])
    boxplot(dat_to_m)
    dat_to_m_t = normalizeQuantiles(dat_to_m)
    boxplot(dat_to_m_t)
    gene_list_df[,7:(dim(gene_list_df)[2]-1)] = dat_to_m_t
  }
  
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
save(comps, file="tmp/comps.RData") # savepoint


# GENE ENRICHMENT ---------------------------------------------------------

load("tmp/total_diff.RData")
load("tmp/total_signal.RData")
load("tmp/total_gene_orders.RData")
load("tmp/fc_res_long_filt.RData")

table(fc_res_long_filt$trmt)
fc_res_long_filt$dir = ifelse(fc_res_long_filt$fc > 0, "up", "down")
fc_res_long_filt$trmt_dir = paste(fc_res_long_filt$trmt, fc_res_long_filt$dir, sep="_")

# up and down separate
# gene_sets = sapply(unique(fc_res_long_filt$trmt_dir), function(x) fc_res_long_filt$gene[fc_res_long_filt$trmt_dir==x])

# up and down together
gene_sets = lapply(unique(fc_res_long_filt$trmt), function(x) unlist(fc_res_long_filt[fc_res_long_filt$trmt==x,'gene']))
names(gene_sets) = unique(fc_res_long_filt$trmt)

lapply(gene_sets, length)

# enrichment

dbs = tbl_df(listEnrichrDbs())
dbs = c("GO_Biological_Process_2017")

enrich_process <- function(x) {
  res = enrichr(unique(x), dbs)[[1]]
  res = tbl_df(filter(res, `Adjusted.P.value` <= 0.05))
}

enrich_res <- lapply(gene_sets, enrich_process)
enrich_res_df = bind_rows(enrich_res, .id = "Trmt")
save(enrich_res, file="tmp/enrich_res.RData")
write_csv(enrich_res_df, "tmp/enrich_res_df.csv")


# MOTIF SEARCH ------------------------------------------------------------

require(RcisTarget)

gene_sets_merged = sapply(unique(fc_res_long_filt$trmt), function(x) fc_res_long_filt$gene[fc_res_long_filt$trmt==x])
lapply(gene_sets_merged, length)

# import motif-gene rankings
# how is this ranking performed?
# how much do the results change varying the tss window?
motif_rankings <- importRankings("~/Downloads/hg19-tss-centered-5kb-7species.mc9nr.feather")

# import motif-tf annotations
data(motifAnnotations_hgnc)

# run enrichment
# question: what motifs are enriched in each cluster?
motif_enrichment <- cisTarget(gene_sets_merged, motif_rankings, motifAnnot=motifAnnotations_hgnc)
table(motif_enrichment$geneSet)

tf_networks = vector("list", length(gene_sets_merged))
names(tf_networks) = names(gene_sets_merged)
tf_df = data.frame()

for(j in 1:length(gene_sets_merged)) { # visualise

  print(j)  
  
  # get the significant motifs
  # can also apply motif-tf  mapping threshold here
  sig_motifs = filter(motif_enrichment, geneSet==names(gene_sets_merged)[j], NES>=4) %>% select(motif) %>% unlist()
  
  # get all interactions between significant tfs and regulated genes
  inc_mat <- getSignificantGenes(gene_sets_merged[[j]], motif_rankings, signifRankingNames=sig_motifs, plotCurve=TRUE, maxRank=5000-20, genesFormat="incidMatrix", method="aprox")$incidMatrix
  motif_mapping = motifAnnotations_hgnc$TF[match(sig_motifs, motifAnnotations_hgnc$motif)]
  rownames(inc_mat) = motif_mapping
  inc_mat = inc_mat[!is.na(rownames(inc_mat)),]
  
  tf_networks[[j]] = list(net=inc_mat, tfs=unique(motif_mapping[!is.na(motif_mapping)]))
  
  edges <- melt(inc_mat)
  edges <- edges[which(edges[,3]==1),1:2]
  edges = distinct(edges)
  names(edges) = c("TF","Gene")
  edges$trmt = names(gene_sets_merged)[j]
  tf_df = rbind(tf_df, edges)
}

# export results
tf_df = tbl_df(left_join(tf_df, fc_res_long_filt, by=c("Gene"="gene","trmt"="trmt")))
write_csv(tf_df, "~/Dropbox/OTAR020/figures_dat/tf_dat.csv")

# upset table
tfs = lapply(tf_networks, function(x) x[[2]])
tfs_table = as.data.frame.matrix((table(stack(tfs))))
tfs_table = cbind(rownames(tfs_table), tfs_table)
rownames(tfs_table) = NULL
upset(tfs_table, order.by="freq", sets=c("primary_macrophage","primary_macrophage_inflamm","THP_1_PMA_LPS","U937_PMA_LPS"))

# construct network 
inc_mat = tf_networks[[which(names(tf_networks)=="THP_1_PMA_LPS")]]$net
edges <- melt(inc_mat)
edges <- edges[which(edges[,3]==1),1:2]
edges = distinct(edges)
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes), label=c(motifs, genes), title=c(motifs, genes), shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))), color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
nodes = distinct(nodes)
nodes = nodes[!duplicated(nodes$id),]

# visualise
visNetwork(nodes, edges) %>% visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE)


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


