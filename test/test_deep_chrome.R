
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
system(paste("ls", paste(c(my_bams_1, my_bams_2)[to_idx], collapse=" "), "| parallel samtools index '{}'"))
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


load("tmp/dc_example.RData")

# add rna expression
dc_example$rna = 0

# genes that are significantly upregulated: pma+lps vs. pma
thp_diff = read_excel("z:/links/bix-analysis-stv/2016/CTTV/THP1/documents/AllGenes_THP1_U937_Altius_DE.xls")
thp_diff$THP1_LPS_log2FoldChange = as.numeric(thp_diff$THP1_LPS_log2FoldChange) 
thp_diff_filt = filter(thp_diff, THP1_LPS_padj<0.05, THP1_LPS_log2FoldChange>1) %>% arrange(desc(THP1_LPS_log2FoldChange))

dc_example$rna[dc_example$Gene %in% thp_diff_filt$gene] = 1
dc_example %>% group_by(Gene) %>% summarise(Status=mean(rna)) %>% dplyr::select(Status) %>% table()

dc_example_norm = dc_example
dc_example_norm[,5:9] = normalizeQuantiles(dc_example[,c(5:9)])

dc_data = matrix(NA, nrow=dim(dc_example_norm)[1]/100, ncol=501)

c_ix = 1
for(i in 1:dim(dc_data)[1]) {
  dc_data[i,1:500] = unlist(dc_example_norm[c_ix:(c_ix+99),c(5:9)])
  dc_data[i,501] = dc_example_norm$rna[c_ix]
  c_ix = c_ix+100
}

dc_data = tbl_df(dc_data)
names(dc_data)[1:500] = as.character(sapply(names(dc_example)[5:9], function(x) paste(x, 1:100, sep="_")))
names(dc_data)[501] = "Y"
dc_data$Y = factor(dc_data$Y)

pheatmap(matrix(as.numeric(dc_data[sample(which(dc_data$Y==1),1),-dim(dc_data)[2]]), nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

dc_data_sample = rbind(
  dc_data[dc_data$Y==1,],
  dc_data[sample(which(dc_data$Y==0),500, replace=FALSE),]
)

task = makeClassifTask(data=dc_data_sample, target="Y")

# lrn = makeLearner("classif.lda")
lrn = makeLearner("classif.randomForest", predict.type="prob", fix.factors.prediction=TRUE)

train_set = sort(c(sample(which(dc_data_sample$Y==1),100, replace=FALSE), sample(which(dc_data_sample$Y==0),250,replace=FALSE)))
test_set = c(1:dim(dc_data_sample)[1])[-train_set]

model = train(lrn, task, subset=train_set)
pred = predict(model, task=task, subset=test_set)
performance(pred, measures=list(mmce, acc))

fv = generateFilterValuesData(task, method="information.gain")
# fv = generateFeatureImportanceData(task, learner=lrn)

pheatmap(matrix(fv$data$information.gain, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE)

av_pos = unlist(apply(filter(dc_data_sample, Y==1) %>% dplyr::select(-Y), 2, mean))
av_neg = unlist(apply(filter(dc_data_sample, Y==0) %>% dplyr::select(-Y), 2, mean))
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

