
load("data/gene_list_all.RData")


# EXPORT COORDINATES ------------------------------------------------------

tss_win = gene_list_all
start(tss_win) = tss_win$transcription_start_site - 5e3
end(tss_win) = tss_win$transcription_start_site + (5e3+1)
write_tsv(data.frame(tss_win)[,1:3], "y:/sandbox/tss_win.bed", col_names=FALSE)


# IMPORT AND MUNGE DATA ---------------------------------------------------

# y = read_tsv("y:/sandbox/bin_regions/otar_samples/thp_1_BR1.txt", col_names=FALSE)
# thp1_df = lapply(as.list(thp1), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/", x), col_names=FALSE))
# y_2 = read_tsv("y:/sandbox/bin_regions/otar_samples/thp_1_BR1_atac.txt", col_names=FALSE)
# y = cbind(y,y_2$X5)
# names(y) = c("Seq","Start","End","Bin","H34K27ac","H3K4me3","CTCF","H3K27me3","ATAC")
# y$Bin = rep(1:100, 18086)
# y$Gene = rep(gene_list_all$hgnc_symbol, each=100)

load("tmp/dc_example.RData")

# add rna expression
dc_example$rna = 0

# genes that are significantly upregulated: pma+lps vs. pma
thp_diff = read_excel("z:/links/bix-analysis-stv/2016/CTTV/THP1/documents/AllGenes_THP1_U937_Altius_DE.xls")
thp_diff$THP1_LPS_log2FoldChange = as.numeric(thp_diff$THP1_LPS_log2FoldChange) 
thp_diff_filt = filter(thp_diff, THP1_LPS_padj<0.05, THP1_LPS_log2FoldChange>1) %>% arrange(desc(THP1_LPS_log2FoldChange))

dc_example$rna[dc_example$Gene %in% thp_diff_filt$gene] = 1
dc_example %>% group_by(Gene) %>% summarise(Status=mean(rna)) %>% dplyr::select(Status) %>% table()

# split into train, test, validate
neg_genes = gene_list_all$hgnc_symbol[!gene_list_all$hgnc_symbol %in% thp_diff_filt$gene]
pos_genes = gene_list_all$hgnc_symbol[gene_list_all$hgnc_symbol %in% thp_diff_filt$gene]

neg_sample = sample(1:3, size=length(neg_genes), replace=TRUE)
pos_sample = sample(1:3, size=length(pos_genes), replace=TRUE)

train_set = rbind(
  dc_example[dc_example$Gene %in% pos_genes[which(pos_sample==1)],],
  dc_example[dc_example$Gene %in% neg_genes[which(neg_sample==1)],]
)

test_set = rbind(
  dc_example[dc_example$Gene %in% pos_genes[which(pos_sample==2)],],
  dc_example[dc_example$Gene %in% neg_genes[which(neg_sample==2)],]
)

validate_set = rbind(
  dc_example[dc_example$Gene %in% pos_genes[which(pos_sample==3)],],
  dc_example[dc_example$Gene %in% neg_genes[which(neg_sample==3)],]
)

# ANALYSIS ----------------------------------------------------------------


