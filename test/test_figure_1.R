
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


# GENERATE DATA -----------------------------------------------------------

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
  dds_list[[i]] <- DESeqDataSet(se_filt[,which(col_data$Mark==marks[i])])
  nrow(dds_list[[i]])
  dds_list[[i]] <- estimateSizeFactors(dds_list[[i]])
  
  # log convert to equalise variance across means (for plotting)  
  rld_list[[i]] = rlog(dds_list[[i]], blind=FALSE)  
  
  # run dds
  dds_list[[i]] <- DESeq(dds_list[[i]])
}


# COMPARE WITH PCA/HEATMAPS -----------------------------------------------

require("pheatmap")
require("RColorBrewer")

for(i in 1:length(rld_list)) {
  
  print(names(rld_list)[i])
  
  # find the equivalent auc data
  auc_ix = which(names(total_data)==names(rld_list[i]))
  my_comp = list(
    count = t(assay(rld_list[[i]])),
    auc = total_data[[auc_ix]]$res[match(
      paste(rld_list[[i]]$Cell, rld_list[[i]]$Stimulus, rld_list[[i]]$Rep, sep ="_"),
      paste(total_data[[auc_ix]]$annot$Cell, total_data[[auc_ix]]$annot$Stimulus, total_data[[auc_ix]]$annot$Rep, sep="_")
    ),]
  )
  
  for(j in 1:length(my_comp)) {
    
    sample_dists <- dist(my_comp[[j]])
    sd_mat <- as.matrix(sample_dists)
    row_names = paste(rld_list[[i]]$Cell, rld_list[[i]]$Stimulus, rld_list[[i]]$Rep, sep ="_")
    rownames(sd_mat) = row_names
    colnames(sd_mat) = NULL
    colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    
    png(paste0("hm_", i, "_", j, ".png"), height=727, width=1052)
    pheatmap(sd_mat, clustering_distance_rows=sample_dists, clustering_distance_cols=sample_dists, col=colors, main=names(rld_list)[i])
    dev.off()
    
    plot_pca(my_comp[[j]], annot_1=row_names, annot_2=rld_list[[i]]$Cell, out_file=paste0("pc_", i, "_", j, ".png"))
    print(plotPCA(rld_list[[i]], intgroup=c("Cell","Stimulus")) + theme_thesis() + ggtitle(names(rld_list)[i]))
    
  }
  
}


# PART 2 - COMPARE COLLAPSE-TO-GENE METHODS -------------------------------

# total_data
# how much difference is there compared to the regulatory genome?

# find dups
dups = duplicated(str_replace(rownames(dat_all[[1]]$H3K27ac$res), "BR[12]_", "")) | duplicated(str_replace(rownames(dat_all[[1]]$H3K27ac$res), "BR[12]_", ""), fromLast=T)
dups = split(which(dups), str_replace(rownames(dat_all[[1]]$H3K27ac$res), "BR[12]_", "")[dups])
marks = names(dat_all[[1]])[1:5]
dist_res = matrix(NA, ncol=length(marks), nrow=3)
colnames(dist_res) = marks
rownames(dist_res) = names(dat_all)[1:3]

for(i in 1:length(marks)) {
  
  for(j in 1:3) { # try each metric
    
    y = dat_all[[j]][[i]]$res
    y_dist = as.matrix(dist(y))
    unlist(lapply(dups, function(x) y_dist[x[1],x[2]]))
    dist_res[j,i] = mean(unlist(lapply(dups, function(x) y_dist[x[1],x[2]])) / mean(y_dist, na.rm=TRUE), na.rm=TRUE)
    
  }
}


