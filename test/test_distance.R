# start off with merged, normalised data - chip + rna
# taken from test_masking.r

# start_data = max reg region over gene body


# LOOP THROUGH RELEVANT LUNG GO TERMS -------------------------------------

require(qusage)
msig_go_bp = read.gmt("c5.all.v6.1.symbols.gmt")

lung_go = c(
  "GO:0006855","GO:0042908","GO:0006281","GO:0004033","GO:0006805","GO:0009494","GO:0017144",
  "GO:0002314","GO:0002377","GO:0002467","GO:0016064","GO:0019882","GO:0042571","GO:0006915",
  "GO:0008219","GO:0000302","GO:0072593","GO:0006802","GO:0006979","GO:0034599","GO:0010193",
  "GO:0070994","GO:0000303","GO:0006801","GO:0006749","GO:0033355","GO:0034635","GO:0036347",
  "GO:1901370","GO:0006809","GO:0071731","GO:0007263","GO:0035713","GO:0046210","GO:0042744",
  "GO:0050665","GO:0045071","GO:0051607","GO:0050691","GO:0045087","GO:0045088","GO:0045089",
  "GO:0002218","GO:0002227","GO:0008063","GO:0034121","GO:0002224","GO:0001816","GO:0001817",
  "GO:0001819","GO:0002775","GO:0002794","GO:0002777","GO:0032602","GO:0032642","GO:0001875",
  "GO:0071222","GO:0038187","GO:1904019"
)

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
go_genes = vector("list", length(lung_go))
names(go_genes) = lung_go

for(i in 1:length(lung_go)) {
  print(paste("Query", i))
  go_res = getBM(attributes=c("hgnc_symbol","go_id"), filters="go", values=lung_go[i], mart=ensembl)
  if(!dim(go_res)[1]) next
  go_genes[[i]] = unique(go_res$hgnc_symbol)
}


# EXAMPLE -----------------------------------------------------------------

load("z:/sandbox/epiChoose/r_data/column_annotation/gene_list_all.RData")

i = 1
genes = go_genes[[i]]
col_ix = which(gene_list_all$hgnc_symbol %in% genes)


  
  end_data = tss_data[[k]]
  
  # slice matrices if necessary
  for(j in 1:length(end_data)) { # each data type
    end_data[[j]]$res = end_data[[j]]$res[,col_ix]
  }
  
  res = dist_mat(end_data, comp_ix=list(c(1:14,18:21), 15), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), my_title=paste(names(go_genes)[i], names(tss_data)[k]), plot_res=FALSE, use_corr=TRUE, font_size=30)
  
  
  
  
  
  
  all_res = matrix(NA, nrow=length(tss_data)*4, ncol=length(go_genes))
  rownames(all_res) = paste(rep(names(tss_data)), rep(c("H3K27ac","H3K4me3","H3K27me3","RNA"), each=6))
  # rownames(all_res) = names(tss_data)
  colnames(all_res) = names(go_genes)
  all_diff = all_res
  r_ix = c(0,6,12,18)
  d_ix = c(1:3,6)
  
  for(i in 1:length(go_genes)) {
    
    print(paste("Building for term", i))
    
    # if(i %in% c(440,573,578,676)) next
    
    # get the genes
    
    
    for(l in 1:4) {
      
      all_diff[k+r_ix[l],i] = mean(unlist(res[d_ix[l],1:2])) - mean(unlist(res[d_ix[l],3:4]))
      if(l==4) all_diff[k+r_ix[l],i] = res[d_ix[l],1] - res[d_ix[l],3]
      if(is.na(all_diff[k+r_ix[l],i])) next
      
      if(all_diff[k+r_ix[l],i] > 0) {
        all_res[k+r_ix[l],i] = "BEAS2B"
      }
      if(all_diff[k+r_ix[l],i] < 0) {
        all_res[k+r_ix[l],i] = "A549"
      }
      
    }
  }
  
}


# PICK OUT BIG DIFFS ------------------------------------------------------

# maximise difference between cell lines
# minimize difference between cell line and cell type

summ_ix = 5
all_res_filt = all_res[summ_ix+c(0,6,12,18),]
all_diff_filt = all_diff[summ_ix+c(0,6,12,18),]
which(abs(all_diff_filt[1,])>0.9 & all_res_filt[1,]=="BEAS2B")


# EXPORT GO TERMS FOR IGV -------------------------------------------------

f <- function(data) {
  nCol <- max(vapply(data, length, 0))
  data <- lapply(data, function(row) c(row, rep(NA, nCol-length(row))))
  data <- matrix(unlist(data), nrow=length(data), ncol=nCol, byrow=TRUE)
  # data = t(data)
  data.frame(data)
}

to_export = f(go_genes)
to_export[is.na(to_export)] = ""
to_export = cbind(names(go_genes), NA, to_export)

write_tsv(to_export, "out.gmt", col_names=FALSE)


# PLOT WITH RNA -----------------------------------------------------------

term_ix = grep("GO_GLYCERALDEHYDE_3_PHOSPHATE_METABOLIC_PROCESS", names(go_genes))

for(term_ix in 29:40) {
  
  if(term_ix %in% c(7,8,9,20,24,28,30,31,36)) next
  
  summ_ix = 5
  
  # get the genes
  genes = go_genes[[term_ix]]
  # genes = c("AKR1C1","BAMBI","ICAM1","KCTD15")
  col_ix = which(gene_list_all$hgnc_symbol %in% genes)
  end_data = tss_data[[summ_ix]]
  
  # slice matrices if necessary
  for(j in 1:length(end_data)) { # each data type
    end_data[[j]]$res = end_data[[j]]$res[,col_ix,drop=FALSE]
  }
  
  png(filename=paste0("final/", str_replace(names(go_genes)[term_ix], ":", "_"), "_plot_1.png"), height=575, width=1271)
  res = dist_mat(end_data, comp_ix=list(c(1:14,18:21), 15), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), my_title=paste("", names(tss_data)[summ_ix]), plot_res=TRUE, use_corr=TRUE, font_size=30)
  dev.off()
  
  x = end_data[[6]]$res[c(1,3,15),]
  x = cbind(x, rownames(x))
  y = melt(x)
  names(y) = c("Cell Line", "Gene", "FPKM")
  png(filename=paste0("final/", str_replace(names(go_genes)[term_ix], ":", "_"), "_plot_2.png"), height=448, width=2014)
  print(ggplot(y, aes(x=Gene, y=FPKM)) + geom_bar(aes(fill=`Cell Line`), position="dodge", stat="identity") + theme_thesis(20) + theme(axis.text.x=element_text(angle=45, hjust=1, size=20)))
  dev.off()
}


# LOOK AT DIFFERENCES -----------------------------------------------------

cut_data = tss_data[[5]] # use max

for(i in 1:length(cut_data)) {
  cut_data[[i]]$res = cut_data[[i]]$res[c(1,3,15),]
  cut_data[[i]]$annot = cut_data[[i]]$annot[c(1,3,15),]
}

spotfire_export = do.call("cbind", lapply(cut_data, function(x) data.frame(log(t(x$res+1)))))
spotfire_export$gene = rownames(spotfire_export)
spotfire_export = spotfire_export[apply(spotfire_export[,grep("RNA", names(spotfire_export))], 1, function(x) !any(is.na(x))),]

# close to a549
spotfire_export$diff_h3k27ac_a549 = abs(spotfire_export$H3K27ac.A549_BR1_Baseline-spotfire_export$H3K27ac.BEAS2B_BR1_Baseline) - abs(spotfire_export$H3K27ac.A549_BR1_Baseline-spotfire_export$H3K27ac.NHBE_BR1_Baseline)
spotfire_export$diff_RNA_a549 = abs(spotfire_export$RNA.A549_BR1_Baseline-spotfire_export$RNA.BEAS2B_BR1_Baseline) - abs(spotfire_export$RNA.A549_BR1_Baseline-spotfire_export$RNA.NHBE_BR1_Baseline)

# close to beas2b
spotfire_export$diff_h3k27ac_beas2b = abs(spotfire_export$H3K27ac.A549_BR1_Baseline-spotfire_export$H3K27ac.BEAS2B_BR1_Baseline) - abs(spotfire_export$H3K27ac.BEAS2B_BR1_Baseline-spotfire_export$H3K27ac.NHBE_BR1_Baseline)
spotfire_export$diff_RNA_beas2b = abs(spotfire_export$RNA.A549_BR1_Baseline-spotfire_export$RNA.BEAS2B_BR1_Baseline) - abs(spotfire_export$RNA.BEAS2B_BR1_Baseline-spotfire_export$RNA.NHBE_BR1_Baseline)

# y$thresh = as.logical(y$BEAS2B_BR1_Baseline<6 & y$A549_BR1_Baseline>9)
# y$thresh[is.na(y$thresh)] = FALSE
# y$label = NA
# y$label[y$thresh] = as.character(y$gene[y$thresh])
# 
# ggplot(y, aes(x=BEAS2B_BR1_Baseline, y=A549_BR1_Baseline)) + geom_point() + theme_thesis() # + geom_text_repel(aes(label=label), fontface="bold")

write_csv(spotfire_export, "tmp/spotfire_export.csv")

# write_tsv(filter(spotfire_export, diff_h3k27ac_a549>1.5 & spotfire_export$diff_RNA_a549>0.9) %>% dplyr::select(gene), "out.txt")


# LOOK AT COURCOT ---------------------------------------------------------

courcot = read_tsv("tmp/courcot_table_1.txt")

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


