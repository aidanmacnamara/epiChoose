
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


# ANALYSIS ----------------------------------------------------------------


