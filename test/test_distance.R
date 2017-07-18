# start off with mask_data

load("r_data/mask_data.RData")
load("r_data/all_data.RData")


# PROJECT 1 RESULT --------------------------------------------------------

start_data = mask_data # choose which dataset

sample_ix = 1:dim(start_data[[1]]$res)[1] # what samples

# sample labels
single_labels = rownames(start_data[[1]]$res)[sample_ix]
group_labels = c(rep("BLUEPRINT",76), rep("GSK",12), rep("ENCODE",18), rep("GSK",31))

# slice matrices if necessary
# for(i in 1:length(start_data)) {
#   start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
#   start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
# }

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=single_labels, marks=names(start_data), plot_type="mds")

png(filename="c:/Downloads/tmp/out.png", height=1200, width=3600)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

res = dist_mat(start_data, comp_ix=list(77:88, 43), labels=single_labels, plot_labels="KU812")


# HOW DO DISTANCES CHANGE WITH GENE SETS? ---------------------------------

start_data = mask_data # choose which dataset

# pick out the lung and blood samples
sample_ix = c(77:88,43,107:113,89)

# sample labels
single_labels = rownames(start_data[[1]]$res)[sample_ix]
group_labels = c(rep("GSK",12), rep("Blueprint",1), rep("GSK",7), rep("ENCODE",1))

# slice matrices if necessary
for(i in 1:length(start_data)) {
  start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
  start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
}

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=single_labels, marks=names(start_data), plot_type="mds")

png(filename="c:/Downloads/tmp/out.png", height=1200, width=3600)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

res = dist_mat(start_data, comp_ix=list(c(1:17,21), 18), labels=single_labels, plot_labels=c("BEAS2B","A549"))


# LOOP THROUGH RELEVANT LUNG GO TERMS -------------------------------------

lung_go = c("GO:0000303","GO:0002314","GO:0002377","GO:0002467","GO:0006749","GO:0006801","GO:0006802","GO:0006809","GO:0006915","GO:0007263","GO:0008219","GO:0010193","GO:0016064","GO:0019882","GO:0033355","GO:0034635","GO:0035713","GO:0036347","GO:0042571","GO:0042744","GO:0046210","GO:0050665","GO:0071731","GO:0072593","GO:1901370")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
term_genes = getBM(attributes=c('hgnc_symbol','go'), filters='go', values=lung_go, mart=ensembl)

load("z:/sandbox/epiChoose/r_data/roi_ensembl_multicell.RData")
load("z:/sandbox/epiChoose/r_data/gene_list_all.RData")

lung_go = names(table(term_genes$go))[table(term_genes$go)>10] # only use those > 10 genes

ranks = c()

for(i in 1:length(lung_go)) {
  
  # get the genes
  genes = filter(term_genes, go==lung_go[i]) %>% dplyr::select(hgnc_symbol) %>% unlist()
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
  res = dist_mat(end_data, comp_ix=list(c(1:17,21), 18), labels=single_labels, plot_labels=c("BEAS2B","A549"), my_title=lung_go[i])
  # dev.off()
  
  ranks = c(
    ranks,
    match(c(18,16,14), order(res [1,])),
    match(c(18,16,14), order(res [2,])),
    match(c(18,16,14), order(res [3,]))
  )
  
}


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
  res = dist_mat(end_data, comp_ix=list(c(1:12,14:21), 13), labels=single_labels, plot_labels=c("F36P","HEL9217","K562","MEG01","UT7","KU812"), my_title=blood_go[i])
  dev.off()
  
  ranks = c(
    ranks,
    match(c(20,15,13), order(res [1,])),
    match(c(20,15,13), order(res [2,])),
    match(c(20,15,13), order(res [3,]))
  )
  
}


# LOOK AT COURCOT ---------------------------------------------------------

courcot = read_tsv("courcot_table_1.txt")

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
  
  png(paste0("c:/Downloads/tmp/", i, ".png"), height=400, width=600)
  res = dist_mat(end_data, comp_ix=list(c(1:12,14:21), 13), labels=single_labels, plot_labels=c("F36P","HEL9217","K562","MEG01","UT7","KU812"), my_title=mono_go[i])
  dev.off()
  
}

