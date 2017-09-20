# start off with mask_data

load("r_data/mask_data_no_bp_with_atac.RData")
load("r_data/mask_data_no_bp_with_atac_col_gene.RData")
load("r_data/all_data.RData")
load("r_data/gsk_chip_filtered.RData")


# TEST WINDOWS ------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC","CTCF") # what marks to select

# get all genes
mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_list_grab = getBM(attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","transcription_start_site"), mart=mart_1, filters=list(chromosome_name=c(as.character(1:22), "X", "Y"), with_protein_id=TRUE))

# choose a tss and tss window
pick_tss <- function(x, p_window=100) {
  
  if(x$strand[1]==1) {
    promoter_end = min(x$transcription_start_site) + p_window
    promoter_start = min(x$transcription_start_site) - p_window
  } else if (x$strand[1]==-1) {
    promoter_start = max(x$transcription_start_site) - p_window
    promoter_end = max(x$transcription_start_site) + p_window
  } else {
    promoter_start = NA
    promoter_end = NA
  }
  
  return(data.frame(cbind(x[1,], promoter_start, promoter_end)))
}

tss_list = list()

for(i in c(500,1e3,2e3)) { # window size around tss
  
  # pick out 1 tss per gene
  gene_list_all = group_by(gene_list_grab, hgnc_symbol) %>% do(pick_tss(., i))
  gene_list_all$strand[gene_list_all$strand==1] = "+"
  gene_list_all$strand[gene_list_all$strand==-1] = "-"
  gene_list_all = makeGRangesFromDataFrame(gene_list_all, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")
  newNames = paste("chr", levels(seqnames(gene_list_all)), sep="")
  names(newNames) = levels(seqnames(gene_list_all))
  gene_list_all = renameSeqlevels(gene_list_all, newNames)
  gene_list_all = gene_list_all[-1]
  gene_list_all = sort(gene_list_all)
  
  my_promoters_gene_all = data.frame(
    seqnames = as.character(seqnames(gene_list_all)),
    start = mcols(gene_list_all)$promoter_start,
    end = mcols(gene_list_all)$promoter_end,
    gene = mcols(gene_list_all)$hgnc_symbol
  )
  my_promoters_gene_all = dplyr::arrange(my_promoters_gene_all, seqnames, start, end) # make sure list is sorted
  tss_list = c(tss_list, list(makeGRangesFromDataFrame(my_promoters_gene_all, keep.extra.columns=TRUE)))
}

# tss_list is a list of granges objects


# GET AUCS ----------------------------------------------------------------

tss_data = vector("list",3)

for(i in 1:length(tss_data)) {
  
  gsk_input = "data/data_gsk.csv"
  gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, tss_list[[i]], marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=4))
  gsk_chip_filtered = prep_across_datatypes(gsk_chip)
  
  encode_input = "data/data_encode.csv"
  encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, tss_list[[i]], marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
  encode_chip_filtered = prep_across_datatypes(encode_chip)
  
  mask_data = vector("list", 5)
  for(i in 1:length(mask_data)) {
    mask_data[[i]]$res = rbind(gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res)
    # renormalize as we are multiple sources
    mask_data[[i]]$res = quantile_norm(mask_data[[i]]$res)
    mask_data[[i]]$annot = bind_rows(gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot)
  }
  
  names(mask_data) = marks
  tss_data[[i]] = mask_data
  
}

# tss_data 3 datasets with different promoter windows


# ADD RNA -----------------------------------------------------------------

p1_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4101-query-results.fpkms.tsv", skip=4)
p3_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4729-query-results.fpkms.tsv", skip=4)

rna = tbl_df(merge(p1_rna, p3_rna, all=TRUE))


# SUMMARISE ENSEMBL REG ---------------------------------------------------

load("r_data/column_annotation/tss_list.RData")
load("r_data/tss_data.RData")

tss_data = c(tss_data, vector("list", 3))

for(j in 1:3) { # merge with rna
  
  start_data = tss_data[[j]] 
  
  # pick out the lung and blood samples
  sample_ix = c(1:17,42:44,72)
  
  # sample labels
  single_labels = rownames(start_data[[1]]$res)[sample_ix]
  group_labels = c(rep("GSK",19), rep("ENCODE",2))
  
  # slice matrices if necessary
  for(i in 1:length(start_data)) {
    start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
    start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
  }
  
  rna_add = prep_rna(rna, tss_list[[j]]$gene, single_labels, c("F36P_BR1_Baseline","HEL9217_BR1_Baseline","K562_BR1_Baseline","KU812_BR1_Baseline","MEG01_BR1_Baseline","UT7_BR1_Baseline","A549_BR1_Baseline","BEAS2B_BR1_Baseline","NHBE_BR1_Baseline"))
  
  start_data[[6]] = rna_add
  names(start_data)[6] = "RNA"
  
  tss_data[[j]] = start_data
  
}

for(j in 4:6) { # collapse reg regions and merge with rna
  
  start_data = mask_data
  
  # pick out the lung and blood samples
  sample_ix = c(1:17,42:44,72)
  
  # sample labels
  single_labels = rownames(start_data[[1]]$res)[sample_ix]
  group_labels = c(rep("GSK",19), rep("ENCODE",2))
  
  # slice matrices if necessary
  for(i in 1:length(start_data)) {
    start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
    start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
  }
  
  rna_add = prep_rna(rna, gene_list_all$hgnc_symbol, single_labels, c("F36P_BR1_Baseline","HEL9217_BR1_Baseline","K562_BR1_Baseline","KU812_BR1_Baseline","MEG01_BR1_Baseline","UT7_BR1_Baseline","A549_BR1_Baseline","BEAS2B_BR1_Baseline","NHBE_BR1_Baseline"))
  
  start_data[[6]] = rna_add
  names(start_data)[6] = "RNA"
  
  # shrink regulatory matrices to gene matrices
  for(i in 1:length(start_data[1:5])) {
    print(paste("Processing data type", names(start_data)[i]))
    if(j==4) {
      start_data[[i]]$res = convert_reg_matrix(start_data[[i]]$res, roi, gene_list_all, reg_window=2000, summ_method="mean")
    }
    if(j==5) {
      start_data[[i]]$res = convert_reg_matrix(start_data[[i]]$res, roi, gene_list_all, reg_window=2000, summ_method="max")
    }
    if(j==6) {
      start_data[[i]]$res = convert_reg_matrix(start_data[[i]]$res, roi, gene_list_all, reg_window=2000, summ_method="sum")
    }
  }
  
  tss_data[[j]] = start_data
  
}

names(tss_data) = c("tss_1","tss_2","tss_4","gene_body_mean","gene_body_max","gene_body_sum")


# ALL DATA SUMMARY --------------------------------------------------------

for(i in 1:length(tss_data)) {
  pca_data = prep_for_plot(tss_data[[i]], annot_1=group_labels, annot_2=single_labels, marks=names(tss_data[[1]]), plot_type="mds")
  
  png(filename=paste0("plots/overall_mds_", names(tss_data)[i], ".png"), height=1200, width=6000)
  print(ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free"))
  dev.off()
  
  png(filename=paste0("plots/overall_dis_", names(tss_data)[i], ".png"), height=575, width=1271)
  res = dist_mat(tss_data[[i]], comp_ix=list(c(1:14,18:20), 15), labels=single_labels, plot_labels=c("BEAS2B","A549"), plot_res=TRUE, use_corr=TRUE)
  dev.off()
}


# LOOP THROUGH RELEVANT LUNG GO TERMS -------------------------------------

lung_go = c(NA,"GO:0006855","GO:0042908","GO:0006281","GO:0004033","GO:0006805","GO:0009494","GO:0017144 ",NA,"GO:0002314","GO:0002377","GO:0002467","GO:0016064","GO:0019882","GO:0042571","GO:0006915","GO:0008219","GO:0000302","GO:0072593","GO:0006802","GO:0006979","GO:0034599","GO:0010193","GO:0070994","GO:0000303","GO:0006801","GO:0006749","GO:0033355","GO:0034635","GO:0036347","GO:1901370","GO:0006809","GO:0071731","GO:0007263","GO:0035713","GO:0046210","GO:0042744","GO:0050665","GO:0045071","GO:0051607")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
term_genes = getBM(attributes=c("hgnc_symbol","go_id"), filters="go", values=lung_go, mart=ensembl)
term_genes = filter(term_genes, go_id %in% lung_go)

# load("z:/sandbox/epiChoose/r_data/column_annotation/roi_ensembl_multicell.RData")
# load("z:/sandbox/epiChoose/r_data/column_annotation/gene_list_all.RData")

# lung_go = names(table(term_genes$go))[table(term_genes$go)>10] # only use those > 10 genes

all_res = matrix(NA, nrow=length(tss_data)*4, ncol=length(lung_go))
rownames(all_res) = paste(rep(names(tss_data)), rep(c("H3K27ac","H3K4me3","H3K27me3","RNA"), each=6))
# rownames(all_res) = names(tss_data)
colnames(all_res) = lung_go
all_diff = all_res

for(i in 1:length(lung_go)) {
  
  if(i %in% c(15)) next
  
  # get the genes
  genes = filter(term_genes, go_id==lung_go[i]) %>% dplyr::select(hgnc_symbol) %>% unlist()
  if(is_empty(genes)) next
  
  # genes_loc = gene_list_all[gene_list_all$hgnc_symbol %in% genes]
  # col_ix = subjectHits(findOverlaps(genes_loc, roi))
  col_ix = which(gene_list_all$hgnc_symbol %in% genes)
  if(length(col_ix)<2) next
  
  for(k in 1:length(tss_data)) { # go through each region list
    
    end_data = tss_data[[k]]
    
    # slice matrices if necessary
    for(j in 1:length(end_data)) { # each data type
      end_data[[j]]$res = end_data[[j]]$res[,col_ix]
    }
    
    res = dist_mat(end_data, comp_ix=list(c(1:14,18:21), 15), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), my_title=paste(lung_go[i], names(tss_data)[k]), plot_res=FALSE, use_corr=TRUE, font_size=30)
    
    r_ix = c(0,6,12,18)
    d_ix = c(1:3,6)
    
    for(l in 1:4) {
      
      all_diff[k+r_ix[l],i] = mean(unlist(res[d_ix[l],1:2])) - mean(unlist(res[d_ix[l],3:4]))
      if(l==4) all_diff[k+r_ix[l],i] = res[d_ix[l],1] - res[d_ix[l],3]
      
      if(all_diff[k+r_ix[l],i] > 0) {
        all_res[k+r_ix[l],i] = "BEAS2B"
      }
      if(all_diff[k+r_ix[l],i] < 0) {
        all_res[k+r_ix[l],i] = "A549"
      }
      
    }
  }
  
}


# EXPORT GO TERMS FOR IGV -------------------------------------------------


u_ids = unique(term_genes$go_id)
to_export = data.frame()

y = tbl_df(term_genes %>% mutate(ind=row_number()) %>% spread(go_id, hgnc_symbol) %>% dplyr::select(-ind))
y[is.na(y)] = ""
y = apply(y, 2, sort)


# LOOK AT DIFFERENCES -----------------------------------------------------

cut_data = start_data

for(i in 1:length(cut_data)) {
  cut_data[[i]]$res = cut_data[[i]]$res[c(1,3,15),]
  cut_data[[i]]$annot = cut_data[[i]]$annot[c(1,3,15),]
}


j=1
y = tbl_df(data.frame(log(t(cut_data[[j]]$res))))
y$gene = rownames(y)

y$thresh = as.logical(y$BEAS2B_BR1_Baseline>10 & y$A549_BR1_Baseline<3)
y$thresh[is.na(y$thresh)] = FALSE
y$label = NA
y$label[y$thresh] = as.character(y$gene[y$thresh])

ggplot(y, aes(x=BEAS2B_BR1_Baseline, y=A549_BR1_Baseline)) + geom_point() + theme_thesis() + geom_text_repel(aes(label=label), fontface="bold")


# LOOK AT COURCOT ---------------------------------------------------------

courcot = read_tsv("c:/Downloads/tmp/courcot_table_1.txt")

genes_loc = gene_list_all[gene_list_all$hgnc_symbol %in% courcot$Gene]

col_ix = subjectHits(findOverlaps(genes_loc, roi))
if(length(col_ix)<2) next

end_data = start_data

# slice matrices if necessary
for(j in 1:length(start_data)) {
  end_data[[j]]$res = start_data[[j]]$res[,col_ix]
}

png(paste0("c:/Downloads/tmp/", i, ".png"), height=400, width=600)
res = dist_mat(end_data, comp_ix=list(c(1:17,21), 18), labels=single_labels, plot_labels=c("BEAS2B","A549"), my_title="Courcot")
dev.off()

# does this correspond to expression differences?
gene_summ = data.frame()
ol = data.frame(findOverlaps(genes_loc, roi))
for(i in 1:length(genes_loc)) {
  gene_summ = rbind(gene_summ, data.frame(Score=mean(start_data[[1]]$res[14,ol[ol$queryHits==i, 'subjectHits']]), Gene = genes_loc$hgnc_symbol[i]))
}

qplot(courcot$`BEAS-2B`, log(gene_summ$Score[match(courcot$Gene, gene_summ$Gene)])) + theme_thesis() + xlab("Expression") + ylab("H3K27ac")


