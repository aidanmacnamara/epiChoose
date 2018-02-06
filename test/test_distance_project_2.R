
require(GenomicAlignments)
require(BiocParallel)
load("data/roi_reg.RData")

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")

# look at counts
# counts of all data types across promoters?
# plan - look at counts first and then try and repeat with auc

files = grep("THP1", list.files("/GWD/bioinfo/projects/RD-Epigenetics-NetworkData/otar_020/GSK/chip-seq/project_2/bam/", full.names=TRUE), value=TRUE)
bamfiles <- BamFileList(files, yieldSize=2000000)
lapply(bamfiles, seqinfo)

register(MulticoreParam())
se <- summarizeOverlaps(features=roi_reg, reads=bamfiles, mode="Union", ignore.strand=TRUE)

col_data = data_gsk[match(colnames(se), str_extract(data_gsk$Bam, "[[:alnum:]\\._-]+$")),] %>% select(Cell, Mark, Rep, Stimulus)
col_data = data.frame(lapply(col_data, factor))

rownames(col_data) = rownames(colData(se))
colData(se) <- DataFrame(col_data)

# slice se to relevant conditions
se_filt = se[,-grep("VD3", colData(se)$Stimulus)]

# some pca

# pma vs. baseline, lps vs. baseline

dds <- DESeqDataSet(se_filt, design=~Mark+Rep+Stimulus)
nrow(dds)
dds <- estimateSizeFactors(dds)

rld = rlog(dds, blind=FALSE)





# RNA RESPONSE ------------------------------------------------------------

thp_diff = read_excel("z:/links/bix-analysis-stv/2016/CTTV/THP1/documents/AllGenes_THP1_U937_Altius_DE.xls")
thp_diff_filt = filter(thp_diff, THP1_LPS_padj<0.01, THP1_LPS_log2FoldChange>2) %>% arrange(desc(THP1_LPS_log2FoldChange))


# PLOTS -------------------------------------------------------------------

n_peaks = unlist(lapply(peak_list, length))
col_data = data.frame(col_data, peaks=n_peaks[sapply(paste(col_data$donor, col_data$time, sep="_"), function(x) grep(x, names(peak_list)))])
ggplot(col_data, aes(x=time, y=peaks)) + geom_boxplot() + theme_thesis()

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

par(mfrow=c(1,2))
plot(log2(counts(dds, normalized=TRUE)[,1:2]+1), pch=16, cex=0.3)
plot(assay(rld)[,1:2], pch=16, cex=0.3)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$cell
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
plotPCA(rld, intgroup=c("time"))


# DESEQ -------------------------------------------------------------------

dds <- DESeq(dds)
res_all = list()
res_all[[1]] <- results(dds, contrast=c("time","3hr","Un"), alpha=0.05, lfcThreshold=1.5)
res_all[[2]] <- results(dds, contrast=c("time","24hr","Un"), alpha=0.05, lfcThreshold=1.5)
res_all[[3]] <- results(dds, contrast=c("time","48hr","Un"), alpha=0.05, lfcThreshold=1.5)


# ANNOTATE ----------------------------------------------------------------

annot_c = annotatePeak(consensus_reduce, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
peak_anno_all = list()
annot_filts = list()

for(i in 1:length(res_all)) {
  
  res = res_all[[i]]
  summary(res)
  
  res = tbl_df(cbind(as.data.frame(annot_c), res))
  p_ix = which(res$padj<=0.05 & res$log2FoldChange>0)
  res_filt = tbl_df(as.data.frame(res)[p_ix,])
  res_filt = res_filt[order(res_filt$log2FoldChange, decreasing=TRUE),]
  
  annot_filt = annotatePeak(consensus_reduce[p_ix], TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
  print(upsetplot(annot_filt))
  
  annot_filts[[i]] = annot_filt
  peak_anno_all[[i]] = res_filt
  
}

venn_func <- function(x) {
  y = paste(as.character(unlist(x[1:3])), collapse="_")
  return(y)
}

require(VennDiagram)
venn.diagram(list(
  hr_3 = apply(peak_anno_all[[1]], 1, venn_func),
  hr_24 = apply(peak_anno_all[[2]], 1, venn_func),
  hr_48 = apply(peak_anno_all[[3]], 1, venn_func)
), filename = "out.tiff")



# ALL DATA SUMMARY --------------------------------------------------------

df_1 = tbl_df(data.frame(Gene=colnames(start_data[[1]]$res), do.call("cbind", lapply(start_data[c(1:3)], function(x) t(x$res)))))
df_rep_1 = dplyr::select(df_1, Gene, contains("BR1")) %>% mutate(Rep=1)
df_rep_2 = dplyr::select(df_1, Gene, contains("BR2")) %>% mutate(Rep=2)

names(df_rep_1) = str_replace_all(names(df_rep_1), "_BR[12]", "")
names(df_rep_2) = str_replace_all(names(df_rep_2), "_BR[12]", "")

df_1 = tbl_df(rbind(df_rep_1, df_rep_2))

names(df_1) = str_replace(names(df_1), "\\.1$", "_H3K4me3")
names(df_1) = str_replace(names(df_1), "\\.2$", "_H3K27me3")
names(df_1) = str_replace(names(df_1), "\\.3$", "_RNA")

df_1_summ = df_1 %>% group_by(Gene) %>% summarise(PMA_H3K27ac=mean(THP.1_PMA), PMA_H3K4me3=mean(THP.1_PMA_H3K4me3), PMA_H3K27me3=mean(THP.1_PMA_H3K27me3))

df_1_summ[,-1] = log(df_1_summ[,-1]+1)
df_1_summ = tbl_df(merge(df_1_summ, thp_res_df, by.x="Gene", by.y="Gene"))
df_1_summ = arrange(df_1_summ, desc(log2FoldChange.x))

# LPS GENE SET ------------------------------------------------------------

load("../epiView/data/msig_go_bp.RData")
m_up = read_tsv("c:/Downloads/tmp/otar_020_tmp/GSE3982_CTRL_VS_LPS_4H_MAC_UP.txt", skip=2, col_names=FALSE)
m_down = read_tsv("c:/Downloads/tmp/otar_020_tmp/GSE3982_CTRL_VS_LPS_4H_MAC_DN.txt", skip=2, col_names=FALSE)

gene_ix = which(df_1_summ$Gene %in% unique(c(m_up$X1, m_down$X1, unlist(msig_go_bp[grep("LIPOPOLYSACCHARIDE", names(msig_go_bp))]))))
# gene_ix = which(df_1$Gene %in% unlist(m_down$X1))

df_1_filtered = df_1_summ # [gene_ix,]

# is there any difference in pma signature between up/non-regulated genes?

df_1_filtered$rna_diff = ifelse(df_1_filtered$log2FoldChange > 0 & df_1_filtered$padj <= 0.01, "Yes", "No")
ggplot(df_1_filtered, aes(x=rna_diff, y=PMA_H3K27ac)) + geom_boxplot() + theme_thesis(20)


# COMBINATIONS ------------------------------------------------------------

ggplot(df_1_filtered, aes(x=PMA_H3K27ac, y=log2FoldChange)) + geom_point(na.rm=TRUE) + theme_thesis(20) + xlab("PMA H3K27ac") + ylab("Expression Change") 

ggplot(df_1_filtered, aes(x=PMA_H3K27ac/PMA_H3K27me3, y=log2FoldChange)) + geom_point(na.rm=TRUE) + theme_thesis(20) + xlab("PMA H3K27ac / PMA H3K27me3") + ylab("Expression Change")

p_1 <- ggplot(df_1_filtered, aes(x=PMA_H3K4me3/PMA_H3K27me3, y=log2FoldChange, text=Gene)) + geom_point(na.rm=TRUE) + theme_thesis(20) + xlab("PMA H3K4me3 / PMA H3K27me3") + ylab("Expression Change")

ggplot(df_1_filtered, aes(x=PMA_H3K27me3/PMA_H3K4me3, y=log2FoldChange)) + geom_point(na.rm=TRUE) + theme_thesis(20) + xlab("PMA H3K27me3 / PMA H3K4me3") + ylab("Expression Change")

ggplotly(p_1, tooltip="Gene")

arrange(df_1_filtered, desc(log2FoldChange))


