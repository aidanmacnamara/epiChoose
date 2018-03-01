
require(GenomicAlignments)
require(BiocParallel)
load("data/roi_reg.RData")

# dat_all - 4 metrics:
# max = largest peak across gene body (gene start/end -+ 2kb)
# tss = largest peak +- 2kb of tss
# sum = sum of peaks across gene body (gene start/end -+ 2kb)
# max 10 closest peaks = not done yet

# METRIC ------------------------------------------------------------------

# test 1: look at quantile normalised aucs versus counts across regulatory regions
# test 2: look at different collapsing-to-gene methods
# metric: test 1; biological reps distance
# metric: test 1/2: biological reps distance and cell type clustering


# TEST 1 ------------------------------------------------------------------

# pick 10*2(?) samples
# for each data type, compare bam vs. auc clustering

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")
files_by_dt = c()
marks = unique(data_gsk$Mark)[1:5]

for(i in 1:length(marks)) {
  mark_filt = filter(data_gsk, file.exists(data_gsk$Bam), Mark==marks[i])
  my_samples = sample(unique(str_replace(mark_filt$Label, "_BR[1-3]", "")), 10, replace=FALSE)
  files_by_dt = c(
    files_by_dt,
    filter(mark_filt, str_replace(Label,  "_BR[1-3]", "") %in% my_samples) %>% dplyr::select(Bam) %>% unlist() %>% as.character() 
  )
}

col_data = data_gsk[match(str_extract(files_by_dt, "[[:alnum:]\\._-]+$"), str_extract(data_gsk$Bam, "[[:alnum:]\\._-]+$")),] %>% dplyr::select(Cell, Mark, Rep, Stimulus)

# get counts over reg regions
bamfiles <- BamFileList(files_by_dt, yieldSize=2000000)
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
