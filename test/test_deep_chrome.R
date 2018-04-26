
load("data/gene_list_all.RData")
data_gsk = read_excel("inst/extdata/data_gsk.xlsx")


# EXPORT COORDINATES ------------------------------------------------------

# this section is making a bed file that splits each gene into 10kb windows with bins of 200bp
tss_win = gene_list_all
start(tss_win) = tss_win$transcription_start_site - 5e3
end(tss_win) = tss_win$transcription_start_site + (5e3+1)
write_tsv(data.frame(tss_win)[,1:3], "ml/tss_win.bed", col_names=FALSE)

# run ml/make_binned_bed.sh


# IMPORT AND MUNGE DATA ---------------------------------------------------

my_label = "THP-1_BR1_Baseline"
my_bams_1 = data_gsk %>% filter(Label==my_label, Mark!="Input", Mark!="ATAC") %>% dplyr::select(Bam) %>% unlist()
my_bams_2 = data_gsk %>% filter(Label==my_label, Mark=="ATAC") %>% dplyr::select(Bam) %>% unlist()

to_idx = !unlist(sapply(c(my_bams_1, my_bams_2), function(x) file.exists(paste(x, ".bai", collapse=" ", sep=""))))
if(any(to_idx)) {
  system(paste("ls", paste(c(my_bams_1, my_bams_2)[to_idx], collapse=" "), "| parallel samtools index '{}'"))
}

my_cmd_1 = paste0("bedtools multicov -bams ", paste(my_bams_1, collapse=" "), " -bed ml/tss_win_binned.bed -f 0.4 > ", "ml/output/tmp_1.txt")
my_cmd_2 = paste0("bedtools multicov -bams ", paste(my_bams_2, collapse=" "), " -bed ml/tss_win_binned.bed -f 0.4 > ", "ml/output/tmp_2.txt")
system(my_cmd_1)
system(my_cmd_2)

# run the code in links/projects/bin_regions

tmp_1 = read_tsv("ml/output/tmp_1.txt", col_names=FALSE)
tmp_2 = read_tsv("ml/output/tmp_2.txt", col_names=FALSE)
tmp_all = cbind(tmp_1, tmp_2$X5)
names(tmp_all) = c("Seq","Start","End","Bin","H34K27ac","H3K4me3","CTCF","H3K27me3","ATAC")
tmp_all$Bin = rep(1:100, 18086)
tmp_all$Gene = rep(gene_list_all$hgnc_symbol, each=100)
assign(my_label, tmp_all)
save(list=my_label, file=paste0("ml/output/", my_label, ".RData"))


# RELOAD DATA -------------------------------------------------------------

load("ml/output/THP-1_BR1_Baseline.RData")
dat = `THP-1_BR1_Baseline`
rm(`THP-1_BR1_Baseline`)


# GET RESPONSE (RNA) ------------------------------------------------------

# run deseq2 first and then filter on non-expressing baseline

rna_dat = read_tsv("inst/extdata/rna/E-MTAB-5191.genes.raw.htseq2.tsv")
rna_annot = read_tsv("inst/extdata/rna/E-MTAB-5191.sdrf.txt")
my_ids = match(names(rna_dat)[-1], rna_annot$`Comment[ENA_RUN]`)
col_data = rna_annot[my_ids,]
col_data$label = col_data$`Source Name`
col_data = separate(col_data, `Source Name`, c("Cell","Rep","Condition"), sep=" ")
col_data$cell_cond = paste(col_data$Cell, col_data$Condition, sep="_")

dds = DESeqDataSetFromMatrix(countData=rna_dat[,-1], colData=col_data, design=~cell_cond)
dds = DESeq(dds)


# PLOT --------------------------------------------------------------------

rld = rlogTransformation(dds)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$label
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
rld$cmpd = str_replace(rld$group, "_.+$", "")
plotPCA(rld, intgroup=c("label"))


# GET DIFFERENTIAL LIST ---------------------------------------------------

# map ensembl ids to symbols
mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), mart=mart_1)
rownames_symbol = mapping$hgnc_symbol[match(rna_dat$`Gene ID`, mapping$ensembl_gene_id)]

res = results(dds, contrast=c('cell_cond','THP1_PMA_RNA','THP1_CTR_RNA'), alpha=0.01)
res = tbl_df(res)
res$ensembl = rna_dat$`Gene ID`
res$symbol = rownames_symbol
res_filt = res %>% filter(padj<1e-8, log2FoldChange>0) %>% arrange(desc(log2FoldChange))

# get genes that have <1 fpkms in control reps
rna_dat_fpkm = read_tsv("inst/extdata/rna/E-MTAB-5191.genes.fpkm.htseq2.tsv")
my_samples = c("THP1 BR1 CTR_RNA","THP1 BR2 CTR_RNA","THP1 BR1 PMA_RNA","THP1 BR2 PMA_RNA")
my_ids = rna_annot$`Comment[ENA_RUN]`[match(my_samples, rna_annot$`Source Name`)]

rna_dat_fpkm_filt = rna_dat_fpkm[,match(my_ids, names(rna_dat_fpkm))]
names(rna_dat_fpkm_filt) = my_samples
rna_dat_fpkm_filt = data.frame(Gene=rna_dat_fpkm$`Gene ID`, rna_dat_fpkm_filt)
# rna_dat_fpkm_filt_long = gather(rna_dat_fpkm_filt, "Sample", "FPKM", 2:(length(my_samples)+1))
# ggplot(rna_dat_fpkm_filt_long, aes(x=FPKM, fill=Sample)) + geom_histogram(bins=2e4) + theme_thesis() + ylab("") + coord_cartesian(xlim=c(0,20))
# table(rna_dat_fpkm_filt_long$FPKM<1)

# what are the genes that are at 0 fpkm pre-stimustim_genes_summ = data.frame(gene=stim_genes$Gene, diff=apply(stim_genes[,4:5]-stim_genes[,2:3], 1, mean, na.rm=TRUE))

rna_dat_fpkm_filt_baseline = rna_dat_fpkm_filt[(rna_dat_fpkm_filt$THP1.BR1.CTR_RNA<1 & rna_dat_fpkm_filt$THP1.BR2.CTR_RNA<1),]
res_filt = res_filt[res_filt$ensembl %in% rna_dat_fpkm_filt_baseline$Gene,]
plot(
  apply(rna_dat_fpkm_filt[match(res_filt$ensembl, rna_dat_fpkm_filt$Gene),2:3], 1, mean),
  apply(rna_dat_fpkm_filt[match(res_filt$ensembl, rna_dat_fpkm_filt$Gene),4:5], 1, mean)
)
abline(h=1)
test_gene = rna_dat_fpkm_filt$Gene[rna_dat_fpkm_filt$THP1.BR1.PMA_RNA<1 & rna_dat_fpkm_filt$THP1.BR2.PMA_RNA<1 & rna_dat_fpkm_filt$Gene %in% res_filt$ensembl][1]
filter(res_filt, ensembl==test_gene)
filter(rna_dat_fpkm_filt, Gene==test_gene)
stim_genes = unique(res_filt$symbol)


# ANALYSIS ----------------------------------------------------------------

# add rna expression
dat$rna = 0
dat$rna[dat$Gene %in% stim_genes] = 1
dat %>% group_by(Gene) %>% summarise(Status=mean(rna)) %>% dplyr::select(Status) %>% table()

# normalize?
dat_norm = dat
dat_norm[,5:9] = normalizeQuantiles(dat[,c(5:9)])

dat_trans = matrix(NA, nrow=dim(dat_norm)[1]/100, ncol=501)

c_ix = 1
for(i in 1:dim(dat_trans)[1]) {
  dat_trans[i,1:500] = unlist(dat_norm[c_ix:(c_ix+99),c(5:9)])
  dat_trans[i,501] = dat_norm$rna[c_ix]
  c_ix = c_ix+100
}

dat_trans = tbl_df(dat_trans)
names(dat_trans)[1:500] = as.character(sapply(names(dat)[5:9], function(x) paste(x, 1:100, sep="_")))
names(dat_trans)[501] = "Y"
dat_trans$Y = factor(dat_trans$Y)

pheatmap(matrix(as.numeric(dat_trans[sample(which(dat_trans$Y==1),1),-dim(dat_trans)[2]]), nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

dat_trans_sample = rbind(
  dat_trans[dat_trans$Y==1,],
  dat_trans[sample(which(dat_trans$Y==0),500, replace=FALSE),]
)

task = makeClassifTask(data=dat_trans_sample, target="Y")

# lrn = makeLearner("classif.lda")
lrn = makeLearner("classif.randomForest", predict.type="prob", fix.factors.prediction=TRUE)

train_set = sort(c(sample(which(dat_trans_sample$Y==1),60, replace=FALSE), sample(which(dat_trans_sample$Y==0),250,replace=FALSE)))
test_set = c(1:dim(dat_trans_sample)[1])[-train_set]

model = train(lrn, task, subset=train_set)
pred = predict(model, task=task, subset=test_set)
performance(pred, measures=list(mmce, acc))

fv = generateFilterValuesData(task, method="information.gain")
# fv = generateFeatureImportanceData(task, learner=lrn)

pheatmap(matrix(fv$data$information.gain, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE)

av_pos = unlist(apply(filter(dat_trans_sample, Y==1) %>% dplyr::select(-Y), 2, mean))
av_neg = unlist(apply(filter(dat_trans_sample, Y==0) %>% dplyr::select(-Y), 2, mean))
av_range = range(c(av_pos, av_neg))

pheatmap(matrix(av_pos, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=av_range[1], to=av_range[2], by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

pheatmap(matrix(av_neg, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=av_range[1], to=av_range[2], by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))


# TRY TOY EXAMPLE ---------------------------------------------------------

require(mlr)

dat_training = read_csv("tmp/train.csv", col_names=FALSE)
dat_testing = read_csv("tmp/test.csv", col_names=FALSE)

x_training = matrix(NA, nrow=20, ncol=500)
x_testing = x_training
c_ix = 1
for(i in 1:dim(x_training)[1]) {
  x_training[i,] = unlist(dat_training[c_ix:(c_ix+99),3:7])
  x_testing[i,] = unlist(dat_testing[c_ix:(c_ix+99),3:7])
  c_ix = c_ix+100
}
x_training = tbl_df(x_training)
x_testing = tbl_df(x_testing)
x_training$Y = factor(rep(c(1,0), each=10))
x_testing$Y = factor(rep(c(1,0), each=10))

pheatmap(matrix(as.numeric(x_training[sample(which(x_training$Y==1),1),-dim(x_training)[2]]), nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

pheatmap(matrix(as.numeric(x_training[sample(which(x_training$Y==0),1),-dim(x_training)[2]]), nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

task = makeClassifTask(data=rbind(x_training, x_testing), target="Y")

# lrn = makeLearner("classif.lda")
lrn = makeLearner("classif.randomForest", predict.type="prob", fix.factors.prediction=TRUE)
train_set = 1:20
test_set = 21:40
model = train(lrn, task, subset=train_set)
pred = predict(model, task=task, subset=test_set)
performance(pred, measures=list(mmce, acc))

