
require(GenomicAlignments)
require(BiocParallel)
load("data/roi_reg.RData")

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")

# look at counts
# counts of all data types across promoters?
# plan - look at counts first and then try and repeat with auc
# /GWD/bioinfo/projects

files = grep("THP1",
             c(
               list.files("z:/links/RD-Epigenetics-NetworkData/otar_020/GSK/chip-seq/project_2/bam/", full.names=TRUE),
               list.files("z:/links/RD-Epigenetics-NetworkData/otar_020/GSK/atac-seq/project_2/bam/", full.names=TRUE)
             ), value=TRUE)

col_data = data_gsk[match(str_extract(files, "[[:alnum:]\\._-]+$"), str_extract(data_gsk$Bam, "[[:alnum:]\\._-]+$")),] %>% dplyr::select(Cell, Mark, Rep, Stimulus)

# take out vds
v_ix = grep("VD3", col_data$Stimulus)
col_data = data.frame(lapply(col_data[-v_ix,], factor))
files = files[-v_ix]

# get counts over reg regions
bamfiles <- BamFileList(files, yieldSize=2000000)
lapply(bamfiles, seqinfo)
register(MulticoreParam())
se <- summarizeOverlaps(features=roi_reg, reads=bamfiles, mode="Union", ignore.strand=TRUE)

rownames(col_data) = rownames(colData(se))
colData(se) <- DataFrame(col_data)

# construct object per mark
marks = unique(col_data$Mark)
dds_list = vector("list", length(marks))
names(dds_list) = marks
rld_list = dds_list

for(i in 1:length(dds_list)) {
  dds_list[[i]] <- DESeqDataSet(se_filt[,which(col_data$Mark==marks[i])], design=~Rep+Stimulus)
  nrow(dds_list[[i]])
  dds_list[[i]] <- estimateSizeFactors(dds_list[[i]])
  
  # log convert to equalise variance across means (for plotting)  
  rld_list[[i]] = rlog(dds_list[[i]], blind=FALSE)  
  
  # run dds
  dds_list[[i]] <- DESeq(dds_list[[i]])
}

require("pheatmap")
require("RColorBrewer")

for(i in 1:length(rld_list)) {
  
  sample_dists <- dist(t(assay(rld_list[[i]])))
  sd_mat <- as.matrix(sample_dists)
  rownames(sd_mat) = paste(rld_list[[i]]$Stimulus, rld_list[[i]]$Rep, sep ="_")
  colnames(sd_mat) = NULL
  colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  png(paste0("hm_", i, ".png"), height=727, width=1052)
  pheatmap(sd_mat, clustering_distance_rows=sample_dists, clustering_distance_cols=sample_dists, col=colors, main=names(rld_list)[i])
  dev.off()
  png(paste0("pc_", i, ".png"), height=624, width=1048)
  print(plotPCA(rld_list[[i]], intgroup=c("Stimulus")) + theme_thesis() + ggtitle(names(rld_list)[i]))
  dev.off()
  
}

# get fpkms
fpkm_list = vector("list", length(dds_list))
names(fpkm_list) = names(dds_list)
for(i in 1:length(fpkm_list)) {
  fpkm_list[[i]] = fpkm(dds_list[[i]])
}


# RNA RESPONSE ------------------------------------------------------------

thp_diff = read_excel("z:/links/bix-analysis-stv/2016/CTTV/THP1/documents/AllGenes_THP1_U937_Altius_DE.xls")
thp_diff$THP1_LPS_log2FoldChange = as.numeric(thp_diff$THP1_LPS_log2FoldChange) 
thp_diff_filt = filter(thp_diff, THP1_LPS_padj<0.01, THP1_LPS_log2FoldChange>1.2) %>% arrange(desc(THP1_LPS_log2FoldChange))

thp_rand = thp_diff[!thp_diff$gene %in% thp_diff_filt$gene & !is.na(thp_diff$THP1_LPS_log2FoldChange) & abs(thp_diff$THP1_LPS_log2FoldChange)<0.5,][
  sample
  (
    1:dim(thp_diff[!thp_diff$gene %in% thp_diff_filt$gene & !is.na(thp_diff$THP1_LPS_log2FoldChange) & abs(thp_diff$THP1_LPS_log2FoldChange)<0.5,])[1],
    dim(thp_diff_filt)[1],
    replace=FALSE
  ),]

gene_model = data.frame(gene=c(thp_diff_filt$gene, thp_rand$gene), rna_lps=c(thp_diff_filt$THP1_LPS_log2FoldChange, thp_rand$THP1_LPS_log2FoldChange), H3K27ac=NA, H3K4me3=NA, H3K27me3=NA, CTCF=NA, ATAC=NA, H3K27ac_LPS=NA, H3K4me3_LPS=NA, H3K27me3_LPS=NA, CTCF_LPS=NA, ATAC_LPS=NA, Group=factor(c(rep("Significant",146), rep("Random",146))))

for(i in 1:length(gene_model$gene)) {
  
  h3k27ac_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model$gene[i])
  gene_model$H3K27ac[i] = mean(fpkm_list$H3K27ac[h3k27ac_ix,7:8])
  gene_model$H3K27ac_LPS[i] = mean(fpkm_list$H3K27ac[h3k27ac_ix,5:6])
  
  h3k4me3_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model$gene[i])
  gene_model$H3K4me3[i] = mean(fpkm_list$H3K4me3[h3k4me3_ix,7:8])
  gene_model$H3K4me3_LPS[i] = mean(fpkm_list$H3K4me3[h3k4me3_ix,5:6])
  
  h3k27me3_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model$gene[i])
  gene_model$H3K27me3[i] = mean(fpkm_list$H3K27me3[h3k27me3_ix,7:8])
  gene_model$H3K27me3_LPS[i] = mean(fpkm_list$H3K27me3[h3k27me3_ix,5:6])
  
  ctcf_ix = which(roi_reg$feature_type_name=="CTCF Binding Site" & roi_reg$SYMBOL==gene_model$gene[i])
  gene_model$CTCF[i] = mean(fpkm_list$CTCF[ctcf_ix,7:8])
  gene_model$CTCF_LPS[i] = mean(fpkm_list$CTCF[ctcf_ix,5:6])
  
  atac_ix = which(roi_reg$feature_type_name=="Open chromatin" & roi_reg$SYMBOL==gene_model$gene[i])
  gene_model$ATAC[i] = mean(fpkm_list$ATAC[atac_ix,7:8])
  gene_model$ATAC_LPS[i] = mean(fpkm_list$ATAC[atac_ix,5:6])
  
}

gene_model_long = gather(gene_model, "Mark", "Mean FPKM", 3:12)
gene_model_long$`Mean FPKM` = log(gene_model_long$`Mean FPKM`+0.1)
gene_model_long$Mark = factor(gene_model_long$Mark)
ggplot(gene_model_long, aes(Group, `Mean FPKM`)) + geom_boxplot() + theme_thesis(20) + facet_wrap(~Mark)
ggplot(gene_model, aes(x=Group, y=rna_lps)) + geom_boxplot() + theme_thesis()

summary(lm(data=gene_model, formula=rna_lps~H3K27ac+H3K4me3+H3K27me3))
summary(glm(data=gene_model, formula=Group~H3K27ac+H3K4me3+H3K27me3, family=binomial(link='logit')))

# is this result consistent with signal data?
load("data/dat_old.RData")

gene_model_signal = gene_model[,c(1:2,13)]
gene_model_signal[,4:6] = NA
names(gene_model_signal)[4:6] = c("H3K27ac","H3K4me3","H3K27me3")

for(i in 1:length(gene_model_signal$gene)) {
  gene_model_signal$H3K27ac[i] = log(mean(dat_all$tss$H3K27ac$res[c(198,204),which(gene_list_all$hgnc_symbol==gene_model_signal$gene[i])]))
  gene_model_signal$H3K4me3[i] = log(mean(dat_all$tss$H3K4me3$res[c(198,204),which(gene_list_all$hgnc_symbol==gene_model_signal$gene[i])]))
  gene_model_signal$H3K27me3[i] = log(mean(dat_all$max$H3K27me3$res[c(198,204),which(gene_list_all$hgnc_symbol==gene_model_signal$gene[i])]))
}

gene_model_signal_long = gather(gene_model_signal, "Mark", "Mean Signal", 4:6)
ggplot(gene_model_signal_long, aes(Group, `Mean Signal`)) + geom_boxplot() + theme_thesis(20) + facet_wrap(~Mark)

y = data.frame(
  Response = thp_diff$THP1_LPS_log2FoldChange,
  H3K27me3 = 
    log(apply(dat_old$H3K27me3$res[c(198,204),match(thp_diff$gene, gene_list_all$hgnc_symbol)], 2, mean)) / 
    log(apply(dat_old$H3K27ac$res[c(198,204),match(thp_diff$gene, gene_list_all$hgnc_symbol)], 2, mean))
)
ggplot(y, aes(H3K27me3, Response)) + geom_point()

par(mfrow=c(2,2))
plot(gene_model$H3K27ac, gene_model_signal$H3K27ac)
plot(gene_model$H3K4me3, gene_model_signal$H3K4me3)
plot(gene_model$H3K27me3, gene_model_signal$H3K27me3)

# pick 10 genes from each group to export to igv
write_tsv(data.frame(gene=sample(gene_model[gene_model$Group=="Random"&!is.nan(gene_model$H3K27ac),'gene'], 40, replace=FALSE)), col_names=FALSE, "out.txt")


# LOOK AT FULL MODEL AND RATIOS -------------------------------------------

gene_model_full = data.frame(gene=thp_diff$gene, rna_lps=thp_diff$THP1_LPS_log2FoldChange, H3K27ac=NA, H3K4me3=NA, H3K27me3=NA, CTCF=NA, ATAC=NA, H3K27ac_LPS=NA, H3K4me3_LPS=NA, H3K27me3_LPS=NA, CTCF_LPS=NA, ATAC_LPS=NA)
gene_model_full = gene_model_full[!is.na(gene_model_full$rna_lps),]
gene_model_full = gene_model_full[sample(1:length(gene_model_full$gene), 1000, replace=FALSE),] # sample

for(i in 1:length(gene_model_full$gene)) {
  
  print(i)
  
  h3k27ac_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model_full$gene[i])
  gene_model_full$H3K27ac[i] = mean(fpkm_list$H3K27ac[h3k27ac_ix,7:8])
  gene_model_full$H3K27ac_LPS[i] = mean(fpkm_list$H3K27ac[h3k27ac_ix,5:6])
  
  h3k4me3_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model_full$gene[i])
  gene_model_full$H3K4me3[i] = mean(fpkm_list$H3K4me3[h3k4me3_ix,7:8])
  gene_model_full$H3K4me3_LPS[i] = mean(fpkm_list$H3K4me3[h3k4me3_ix,5:6])
  
  h3k27me3_ix = which(grepl("Promoter|Enhancer", roi_reg$feature_type_name) & roi_reg$SYMBOL==gene_model_full$gene[i])
  gene_model_full$H3K27me3[i] = mean(fpkm_list$H3K27me3[h3k27me3_ix,7:8])
  gene_model_full$H3K27me3_LPS[i] = mean(fpkm_list$H3K27me3[h3k27me3_ix,5:6])
  
  ctcf_ix = which(roi_reg$feature_type_name=="CTCF Binding Site" & roi_reg$SYMBOL==gene_model_full$gene[i])
  gene_model_full$CTCF[i] = mean(fpkm_list$CTCF[ctcf_ix,7:8])
  gene_model_full$CTCF_LPS[i] = mean(fpkm_list$CTCF[ctcf_ix,5:6])
  
  atac_ix = which(roi_reg$feature_type_name=="Open chromatin" & roi_reg$SYMBOL==gene_model_full$gene[i])
  gene_model_full$ATAC[i] = mean(fpkm_list$ATAC[atac_ix,7:8])
  gene_model_full$ATAC_LPS[i] = mean(fpkm_list$ATAC[atac_ix,5:6])
  
}

ggplot(gene_model_full, aes(x=H3K27ac, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=H3K27me3, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=H3K27ac/H3K27me3, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=H3K27me3/H3K27ac, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=H3K27ac/H3K4me3, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=H3K27me3/H3K4me3, y=rna_lps)) + geom_point() + theme_thesis()
ggplot(gene_model_full, aes(x=ATAC, y=rna_lps)) + geom_point() + theme_thesis()


# RNA ABSOLUTE VALUES -----------------------------------------------------

thp1 = list.files("~/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/")
thp1_df = lapply(as.list(thp1), function(x) read_tsv(paste0("~/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/", x), col_names=FALSE))
# merge into sample per column
my_names = thp1_df[[1]]$X1
rna = data.frame(do.call("cbind", lapply(thp1_df, function(x) x$X2)))
rna = cbind(my_names, rna)

# get names
sample_info = read_csv("z:/links/bix-analysis-stv/2016/CTTV/U937/data/sampleInfo.csv", col_names=FALSE)
all_names = sample_info$X1[match(str_extract(thp1, "^[[:alnum:]]+"), sample_info$X32)]
names(rna) = c("Gene Name", all_names)



