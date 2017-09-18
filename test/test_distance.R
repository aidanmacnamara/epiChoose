# start off with mask_data

load("r_data/mask_data_no_bp_with_atac.RData")
load("r_data/mask_data_no_bp_with_atac_col_gene.RData")
load("r_data/all_data.RData")
load("r_data/gsk_chip_filtered.RData")


# TEST WINDOWS ------------------------------------------------------------

marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC","CTCF")

# link to ensembl
mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_list_grab = getBM(attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","transcription_start_site"), mart=mart_1, filters=list(chromosome_name=c(as.character(1:22), "X", "Y"), with_protein_id=TRUE))

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


# HOW DO DISTANCES CHANGE WITH GENE SETS? ---------------------------------

start_data = mask_data # choose which dataset

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

### incorporate rna data ###

p1_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4101-query-results.fpkms.tsv", skip=4)
p3_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4729-query-results.fpkms.tsv", skip=4)

rna = tbl_df(merge(p1_rna, p3_rna, all=TRUE))

rna_add = prep_rna(rna, gene_list_all$hgnc_symbol, single_labels, c("F36P_BR1_Baseline","HEL9217_BR1_Baseline","K562_BR1_Baseline","KU812_BR1_Baseline","MEG01_BR1_Baseline","UT7_BR1_Baseline","A549_BR1_Baseline","BEAS2B_BR1_Baseline","NHBE_BR1_Baseline"))

start_data[[6]] = rna_add
names(start_data)[6] = "RNA"

# shrink regulatory matrices to gene matrices
for(i in 1:length(start_data[1:5])) {
  print(paste("Processing data type", names(start_data)[i]))
  start_data[[i]]$res = convert_reg_matrix(start_data[[i]]$res, roi, gene_list_all, reg_window=2000, summ_method="max")
}

plot(log(start_data[[1]]$res[1,]), log(start_data[[6]]$res[1,]))

### /incorporate rna data/ ###

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=single_labels, marks=names(start_data), plot_type="mds")

png(filename="tmp/out.png", height=1200, width=6000)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

res = dist_mat(start_data, comp_ix=list(c(1:14,18:20), 15), labels=single_labels, plot_labels=c("BEAS2B","A549"), plot_res=TRUE, use_corr=TRUE)


# LOOP THROUGH RELEVANT LUNG GO TERMS -------------------------------------

lung_go = c("GO:0000303","GO:0002314","GO:0002377","GO:0002467","GO:0006749","GO:0006801","GO:0006802","GO:0006809","GO:0006915","GO:0007263","GO:0008219","GO:0010193","GO:0016064","GO:0019882","GO:0033355","GO:0034635","GO:0035713","GO:0036347","GO:0042571","GO:0042744","GO:0046210","GO:0050665","GO:0071731","GO:0072593","GO:1901370","GO:0045071","GO:0051607")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
term_genes = getBM(attributes=c("hgnc_symbol","go_id"), filters="go", values=lung_go, mart=ensembl)
term_genes = filter(term_genes, go_id %in% lung_go)

load("z:/sandbox/epiChoose/r_data/column_annotation/roi_ensembl_multicell.RData")
load("z:/sandbox/epiChoose/r_data/column_annotation/gene_list_all.RData")

lung_go = names(table(term_genes$go))[table(term_genes$go)>10] # only use those > 10 genes
r_squared = rep(NA, length(lung_go))

for(i in 1:length(lung_go)) {
  
  if(i==10) next
  
  # get the genes
  genes = filter(term_genes, go_id==lung_go[i]) %>% dplyr::select(hgnc_symbol) %>% unlist()
  if(is_empty(genes)) next
  
  # genes_loc = gene_list_all[gene_list_all$hgnc_symbol %in% genes]
  # col_ix = subjectHits(findOverlaps(genes_loc, roi))
  col_ix = which(gene_list_all$hgnc_symbol %in% genes)
  if(length(col_ix)<2) next
  
  end_data = start_data
  
  # slice matrices if necessary
  for(j in 1:length(start_data)) {
    end_data[[j]]$res = start_data[[j]]$res[,col_ix]
  }
  
  png(paste0("tmp/plot_", i, "_1_max.png"), height=575, width=1271)
  res = dist_mat(end_data, comp_ix=list(c(1:14,18:21), 15), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), my_title=lung_go[i], plot_res=TRUE, use_corr=TRUE, font_size=30)
  dev.off()
  
  x = end_data[[6]]$res[c(1,3,15),]
  x = cbind(x, rownames(x))
  y = melt(x)
  names(y) = c("Cell Line", "Gene", "FPKM")
  png(paste0("tmp/plot_", i, "_2_max.png"), height=448, width=2014)
  print(ggplot(y, aes(x=Gene, y=FPKM)) + geom_bar(aes(fill=`Cell Line`), position="dodge", stat="identity") + theme_thesis(10) + theme(axis.text.x=element_text(angle=45, hjust=1, size=10)))
  dev.off()
  png(paste0("tmp/plot_", i, "_3_max.png"), height=800, width=900)
  print(ggplot(y, aes(x=`Cell Line`, y=log(FPKM))) + geom_boxplot() + theme_thesis() + theme(axis.text.x=element_text(angle=45, hjust=1)))
  dev.off()
  
  lm_input = data.frame(
    rna_response = as.numeric(t(end_data[[6]]$res[c(1,3,15),])),
    k27ac = as.numeric(t(end_data[[1]]$res[c(1,3,15),])),
    k4me3 = as.numeric(t(end_data[[2]]$res[c(1,3,15),])),
    k27me3 = as.numeric(t(end_data[[3]]$res[c(1,3,15),])),
    # atac = as.numeric(t(end_data[[4]]$res[c(1,3,15),])),
    cell_type = factor(rep(start_data[[1]]$annot$Label[c(1,3,15)], each=dim(end_data[[6]]$res[c(1,3,15),])[2])),
    gene = factor(rep(colnames(end_data[[1]]$res), dim(end_data[[1]]$res)[1]))
  )
  
  lm_res = lm(rna_response ~ k27ac + k4me3 + k27me3 + k27ac:k4me3:k27me3, data=lm_input)
  r_squared[i] = summary(lm_res)$adj.r.squared
  
  png(paste0("tmp/plot_", i, "_4_max.png"), height=632, width=1446)
  print(ggpairs(dplyr::select(lm_input, rna_response, k27ac, k4me3, k27me3)) + theme_thesis(20))
  dev.off()
  
}


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


# LOOP THROUGH BLOOD GO TERMS ---------------------------------------------

blood_go = c("GO:0045646","GO:0030218","GO:0060319","GO:0043249","GO:0048821","GO:0043362","GO:0030219","GO:0035855","GO:0045652","GO:0036017","GO:0036018","GO:0038162","GO:0042977","GO:0010534","GO:0010535","GO:1902569","GO:0030097","GO:0035166","GO:0071425","GO:0030099","GO:0002244","GO:0060218","GO:1903706","GO:0051607","GO:0042908","GO:0052697","GO:0006805","GO:0042178","GO:0009410","GO:0042910","GO:0071466","GO:0090484","GO:0046618","GO:0008144","GO:0015893","GO:0042493","GO:0017144","GO:0042737","GO:0035690","GO:0090484","GO:0015695","GO:0007254","GO:0046328","GO:0004705","GO:0038066","GO:0004896","GO:0051607","GO:0030324","GO:0060428")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
term_genes = getBM(attributes=c('hgnc_symbol','go'), filters='go', values=blood_go, mart=ensembl)

blood_go = names(table(term_genes$go))[table(term_genes$go)>10] # only use those > 10 genes

ranks = c()

for(i in 1:length(blood_go)) {
  
  # get the genes
  genes = filter(term_genes, go==blood_go[i]) %>% dplyr::select(hgnc_symbol) %>% unlist()
  if(is_empty(genes)) next
  
  genes_loc = gene_list_all[gene_list_all$hgnc_symbol %in% genes]
  
  col_ix = subjectHits(findOverlaps(genes_loc, roi))
  if(length(col_ix)<2) next
  
  end_data = start_data
  
  # slice matrices if necessary
  for(j in 1:length(start_data)) {
    end_data[[j]]$res = start_data[[j]]$res[,col_ix]
  }
  
  png(paste0("c:/Downloads/tmp/", i, ".png"), height=400, width=600)
  res = dist_mat(end_data, comp_ix=list(c(1:12,14:21), 13), labels=single_labels, plot_labels=c("F36P","HEL9217","K562","MEG01","UT7","KU812"), my_title=blood_go[i], plot_res=FALSE)
  dev.off()
  
  ranks = c(
    ranks,
    match(c(20,15,13), order(res [1,])),
    match(c(20,15,13), order(res [2,])),
    match(c(20,15,13), order(res [3,]))
  )
  
}


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


# MONOCYTE MACROPHAGE MODELS ----------------------------------------------

start_data = mask_data # choose which dataset

# pick out the monocyte/macrophage samples
sample_ix = c(grep("monocyte|macrophage", rownames(mask_data[[1]]$res)), 114:137, 101)

# sample labels
single_labels = rownames(start_data[[1]]$res)[sample_ix]
single_labels[1:20] = str_replace(single_labels[1:20], "[[:alnum:]]+_", "")
group_labels = c(rep("Blueprint",20), rep("GSK",24), rep("ENCODE",1))

# slice matrices if necessary
for(i in 1:length(start_data)) {
  start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
  start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
}

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=single_labels, marks=names(start_data), plot_type="mds")

png(filename="c:/Downloads/tmp/out.png", height=1200, width=3600)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

mono_go = c("GO:0002548","GO:0042117","GO:0035702","GO:0070487","GO:0035696","GO:0061516","GO:0030224")
term_genes = getBM(attributes=c('hgnc_symbol','go'), filters='go', values=mono_go, mart=ensembl)

mono_go = names(table(term_genes$go))[table(term_genes$go)>10] # only use those > 10 genes

for(i in 1:length(mono_go)) {
  
  # get the genes
  genes = filter(term_genes, go==mono_go[i]) %>% dplyr::select(hgnc_symbol) %>% unlist()
  if(is_empty(genes)) next
  
  genes_loc = gene_list_all[gene_list_all$hgnc_symbol %in% genes]
  
  col_ix = subjectHits(findOverlaps(genes_loc, roi))
  if(length(col_ix)<2) next
  
  end_data = start_data
  
  # slice matrices if necessary
  for(j in 1:length(start_data)) {
    end_data[[j]]$res = start_data[[j]]$res[,col_ix]
  }
  
  # png(paste0("c:/Downloads/tmp/", i, ".png"), height=400, width=600)
  res = dist_mat(end_data, comp_ix=list(c(21:44), list(c(1:4,6,45))), labels=single_labels, plot_labels=c("THP-1","U937"), my_title=paste("Monocytes:", mono_go[i]))
  res = dist_mat(end_data, comp_ix=list(c(21:44), list(c(5,7:20))), labels=single_labels, plot_labels=c("THP-1","U937"), my_title=paste("Macrophages:", mono_go[i]))
  # dev.off()
  
}

