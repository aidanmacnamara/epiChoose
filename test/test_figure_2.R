# comparison of monocytes/macrophages with u937/thp-1

# 1. establish the amount of paired samples in blueprint
# 2. groups: monocytes, macrophages, u937, thp-1
# 3. how to define dynamic regions?
# 4. figure 4d for dynamic regions i.e. regions that differ between monocytes and macrophages and what they look like in u937/thp-1
# 5. tf inference from deregulated sites (look at luz's review for methods)


# PAIRED SAMPLES IN BLUEPRINT ---------------------------------------------

# get the necessary files from blueprint
bp <- read_tsv("inst/extdata/blueprint_files.tsv")
bp[is.na(bp)] = "" # removes nas that affect dplyr selections

# which donors have monocyte and macrophage samples?
bp_tbl = bp %>% group_by(Donor) %>% summarise("Both"=(any(grepl("monocyte",`Sub-group`,ignore.case=TRUE))) & (any(grepl("macrophage",`Sub-group`,ignore.case=TRUE))))
bp_filt = filter(bp, Donor %in% unlist(bp_tbl[bp_tbl$Both==TRUE, 'Donor']), Format=="bigWig")
bp_filt %>% group_by(Donor) %>% summarise(Marks=all(c("H3K4me3","H3K27ac") %in% Experiment))


# DYNAMIC H3K27AC REGIONS -------------------------------------------------

load("data/total_data.RData")
load("data/dat_all.RData")

dat = dat_all$tss$H3K27ac$res[grep("monocyte|macrophage", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE),]
rownames(dat) = str_extract(rownames(dat), "[Mm]onocyte|[Mm]acrophage")
rownames(dat) = tolower(rownames(dat))
rownames(dat) = str_replace(rownames(dat), "_BR[12]", "")

# remove na columns
dim(dat)
na_ix = which(apply(is.na(dat), 2, all))
dat = dat[,-na_ix] # remove regions with no data
# dat = dat[-which(apply(dat, 1, function(x) all(is.na(x)))),]  # remove samples with no data
dim(dat)

# pick high-variance regions
head(order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE))
ggplot(data.frame(Type=names(dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), AUC=dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), aes(Type, AUC)) + geom_boxplot() + theme_thesis(20)
col_ix = head(order(apply(dat,2,var,na.rm=TRUE), decreasing=TRUE), 1000)
dat_filt = dat[,col_ix]
pheatmap(log(dat_filt), show_colnames=FALSE)

# kmeans
dat_filt_t = t(dat_filt)
wss <- (nrow(dat_filt_t)-1) * sum(apply(dat_filt_t,2,var))
for(i in 2:15) {
  wss[i] <- sum(kmeans(dat_filt_t, centers=i)$withinss)
}

plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within-group sum-of-squares")
fit <- kmeans(dat_filt_t, 7)
col_clust = hclust(dist(dat_filt_t))

pheatmap(log(dat_filt),
         annotation_col = data.frame(Cluster=factor(fit$cluster)),
         cluster_cols = col_clust,
         show_colnames = FALSE
)

which_clust = 7
to_plot = rbind(
  data.frame(Type="Macrophage", AUC=as.numeric(dat_filt[rownames(dat_filt)=="macrophage",fit$cluster==which_clust])),
  data.frame(Type="Monocyte", AUC=as.numeric(dat_filt[rownames(dat_filt)=="monocyte",fit$cluster==which_clust]))
)
ggplot(to_plot, aes(Type, log(AUC))) + geom_boxplot() + theme_thesis(20)

# add thp-1 and u937
dat_add = dat_all$tss$H3K27ac$res[grep("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE),]
dim(dat_add)
dat_add = dat_add[,-na_ix]
dim(dat_add)[2] == dim(dat)[2]
dat_add_filt = dat_add[,col_ix]
dim(dat_add_filt)
all(colnames(dat_filt) == colnames(dat_add_filt))
dat_add_filt = dat_add_filt[-which(apply(dat_add_filt, 1, function(x) all(is.na(x)))),]
dim(dat_add_filt)

pheatmap(log(rbind(dat_filt, dat_add_filt)),
         annotation_col = data.frame(Cluster=factor(fit$cluster)),
         cluster_cols = col_clust,
         show_colnames = FALSE
)

# what are the clusters?

write_tsv(data.frame(names(which(fit$cluster==7))), "out_1.txt", col_names=FALSE) # monocyte high
# http://amp.pharm.mssm.edu/Enrichr/enrich?dataset=42vf0

write_tsv(data.frame(names(which(fit$cluster==5))), "out_2.txt", col_names=FALSE) # macrophage high
# http://amp.pharm.mssm.edu/Enrichr/enrich?dataset=42vf4


# APPLY T-TEST ------------------------------------------------------------

p_res_g = rep(NA, dim(dat_filt)[2])
for(i in 1:length(p_res_g)) {
  p_res_g[i] = t.test(dat_filt[which(rownames(dat_filt)=="monocyte"),i], dat[which(rownames(dat_filt)=="macrophage"),i], alternative="g")$p.value
}

check_ix = 1
to_plot = rbind(
  data.frame(Type="Macrophage", AUC=as.numeric(dat_filt[rownames(dat_filt)=="macrophage", order(p_res_g)[check_ix]])),
  data.frame(Type="Monocyte", AUC=as.numeric(dat_filt[rownames(dat_filt)=="monocyte", order(p_res_g)[check_ix]]))
)
ggplot(to_plot, aes(Type, log(AUC))) + geom_boxplot() + theme_thesis(20) + ggtitle(colnames(dat_filt)[order(p_res_g)[check_ix]])


# LOG FOLD CHANGE CALCULATION ---------------------------------------------

load("data/dat_all.RData")
load("data/gene_list_all.RData")
tss_window = 2e3
tss_regions = gene_list_all
start(tss_regions) = tss_regions$transcription_start_site - tss_window
end(tss_regions) = tss_regions$transcription_start_site + tss_window


# COUNTS ------------------------------------------------------------------

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")
bamfiles = filter(data_gsk, grepl("thp", Label, ignore.case=TRUE), Mark=="H3K27ac") %>% dplyr::select(Bam) %>% unlist()

col_data = data_gsk[match(str_extract(bamfiles, "[[:alnum:]\\._-]+$"), str_extract(data_gsk$Bam, "[[:alnum:]\\._-]+$")),] %>% dplyr::select(Cell, Mark, Rep, Stimulus)
col_data = as.data.frame(col_data)

bamfiles <- BamFileList(bamfiles, yieldSize=2000000)
lapply(bamfiles, seqinfo)
register(MulticoreParam())
se <- summarizeOverlaps(features=tss_regions, reads=bamfiles, mode="Union", ignore.strand=TRUE)

rownames(col_data) = rownames(colData(se))
colData(se) <- DataFrame(col_data)

dds = DESeqDataSet(se, design=~Stimulus)
dds = DESeq(dds)
rld = rlog(dds, blind=FALSE) 

res = results(dds, contrast=c("Stimulus","PMA","Baseline"))
res$hgnc_symbol = gene_list_all$hgnc_symbol
res = tbl_df(res)
res_filt = filter(res, padj < 0.05, abs(log2FoldChange) > 1.2) %>% arrange(desc(abs(log2FoldChange)))


# AUC ---------------------------------------------------------------------

# 1. start off with h3k27ac matrix
dat_add = dat_all$tss$H3K27ac$res[grep("thp-1", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE),]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)

col_data_auc = col_data
rownames(col_data_auc) = rownames(dat_add)
dds_auc = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_auc, design=~Stimulus)

dds_auc = DESeq(dds_auc)
rld_auc = rlog(dds_auc, blind=FALSE) 

res_auc = results(dds_auc, contrast=c("Stimulus","PMA","Baseline"))
res_auc$hgnc_symbol = gene_list_all$hgnc_symbol
res_auc = tbl_df(res_auc)
res_auc_filt = filter(res_auc, padj < 0.05, abs(log2FoldChange) > 1.2) %>% arrange(desc(abs(log2FoldChange)))

res_filt$auc_diff = res_auc$log2FoldChange[match(res_filt$hgnc_symbol,res_auc$hgnc_symbol)]
plot(res_filt$log2FoldChange, res_auc_filt$log2FoldChange[match(res_filt$hgnc_symbol, res_auc_filt$hgnc_symbol)])
plot(euler(list(AUC=res_auc_filt$hgnc_symbol, Counts=res_filt$hgnc_symbol)), quantities=TRUE)

all_genes = unique(c(res_filt$hgnc_symbol, res_auc_filt$hgnc_symbol))
plot(
  res_filt$log2FoldChange[match(all_genes, res_filt$hgnc_symbol)],
  res_auc_filt$log2FoldChange[match(all_genes, res_auc_filt$hgnc_symbol)]
)



# order samples for visualisation
s_order = c(1,7,2,8,3,9,4,10,5,11,6,12)
dat_add = dat_add[s_order,]

# 2. establish dynamic regions
ctrl_ix = which(grepl("THP-1_BR[12]_Baseline", rownames(dat_add)))
trmt_ix = which(grepl("THP-1_BR[12]_PMA$", rownames(dat_add)))

res_auc = data.frame(p_value=rep(NA,length(gene_list_all)), fc=rep(NA,length(gene_list_all)), gene=gene_list_all$hgnc_symbol)

get_p_value <- function(x) {
  if(any(is.na(c(x[trmt_ix],x[ctrl_ix]))) | any(c(x[trmt_ix],x[ctrl_ix])==0)) {
    return(NA)
  } else {
    return(t.test(log(x[trmt_ix]+0.1),log(x[ctrl_ix]+0.1))$p.value)
  }
}

get_fc <- function(x) {
  if(any(is.na(c(x[trmt_ix],x[ctrl_ix])))) {
    return(NA)
  } else {
    if(mean(x[trmt_ix]) > mean(x[ctrl_ix])) {
      return(log(mean(x[trmt_ix])+0.1)-log(mean(x[ctrl_ix])+0.1))
    } else {
      return(-(log(mean(x[ctrl_ix])+0.1)-log(mean(x[trmt_ix])+0.1)))
    }
  }
}  

res_auc$p_value = apply(dat_add, 2, get_p_value)
res_auc$fc = apply(dat_add, 2, get_fc)
res_auc %>% arrange(desc(abs(fc))) %>% head

