
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

require(BiocParallel)
marks = c("H3K27ac","H3K4me3","H3K27me3","CTCF")

my_labels = c("THP-1_BR1_Baseline","THP-1_BR2_Baseline","THP-1_BR1_PMA","THP-1_BR2_PMA","U937_BR1_Baseline","U937_BR2_Baseline","U937_BR1_PMA","U937_BR2_PMA")

ml_generate <- function(x) {
  
  my_bams_1 = data_gsk %>% filter(Label==my_labels[x], Mark!="Input", Mark!="ATAC", Mark!="Control") %>% dplyr::select(Bam) %>% unlist()
  
  # make sure files are ordered
  my_bams_1 = my_bams_1[match(marks, sapply(my_bams_1, function(x) str_extract(x, paste(marks, collapse="|"))))]
  my_bams_2 = data_gsk %>% filter(Label==my_labels[x], Mark=="ATAC") %>% dplyr::select(Bam) %>% unlist()
  
  to_idx = !unlist(sapply(c(my_bams_1, my_bams_2), function(x) file.exists(paste(x, ".bai", collapse=" ", sep=""))))
  if(any(to_idx)) {
    system(paste("ls", paste(c(my_bams_1, my_bams_2)[to_idx], collapse=" "), "| parallel samtools index '{}'"))
  }
  
  tmp_1 = str_replace(tempfile(tmpdir=""), "^\\\\", "")
  tmp_2 = str_replace(tempfile(tmpdir=""), "^\\\\", "")
  
  my_cmd_1 = paste0("bedtools multicov -bams ", paste(my_bams_1, collapse=" "), " -bed ml/tss_win_binned.bed -f 0.4 > ", "ml/output/", tmp_1, ".txt")
  my_cmd_2 = paste0("bedtools multicov -bams ", paste(my_bams_2, collapse=" "), " -bed ml/tss_win_binned.bed -f 0.4 > ", "ml/output/", tmp_2, ".txt")
  system(my_cmd_1)
  system(my_cmd_2)
  
  # run the code in links/projects/bin_regions
  
  tmp_1 = read_tsv(paste0("ml/output/", tmp_1, ".txt"), col_names=FALSE)
  tmp_2 = read_tsv(paste0("ml/output/", tmp_2, ".txt"), col_names=FALSE)
  tmp_all = cbind(tmp_1, tmp_2$X5)
  names(tmp_all) = c("Seq","Start","End","Bin",marks,"ATAC")
  tmp_all$Bin = rep(1:100, 18086)
  tmp_all$Gene = rep(gene_list_all$hgnc_symbol, each=100)
  assign(my_labels[x], tmp_all)
  save(list=my_labels[x], file=paste0("ml/output/", my_labels[x], ".RData"))
  system(paste0("rm ml/output/", tmp_1, ".txt"))
  system(paste0("rm ml/output/", tmp_2, ".txt"))
}

bplapply(seq(along=my_labels), ml_generate, BPPARAM=MulticoreParam(workers=length(my_labels)))


# RELOAD DATA -------------------------------------------------------------

# my_files = list.files("ml/output/")
# my_data = vector("list", length(my_files))
# names(my_data) = str_replace(my_files, "\\.RData", "")
# for(i in 1:length(my_data)) {
#   my_data[[i]] = get(load(paste0("ml/output/", my_files[[i]])))
# }
load("ml/output/my_data.RData")


# GET RESPONSE (RNA) ------------------------------------------------------

# run deseq2 first and then filter on non-expressing baseline

rna_dat = read_tsv("inst/extdata/rna/E-MTAB-5191.genes.raw.htseq2.tsv") # get data
rna_annot = read_tsv("inst/extdata/rna/E-MTAB-5191.sdrf.txt") # get annotation
my_ids = match(names(rna_dat)[-1], rna_annot$`Comment[ENA_RUN]`) # match data/annotation ids

# construct sample table
col_data = rna_annot[my_ids,]
col_data$label = col_data$`Source Name`
col_data = separate(col_data, `Source Name`, c("Cell","Rep","Condition"), sep=" ")
col_data$cell_cond = paste(col_data$Cell, col_data$Condition, sep="_")

# run deseq2
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

res_diff = vector("list", 4)
names(res_diff) = c("THP-1_Baseline","THP-1_PMA","U937_Baseline","U937_PMA")
res_diff[[1]] = results(dds, contrast=c('cell_cond','THP1_PMA_RNA','THP1_CTR_RNA'), alpha=0.01)
res_diff[[2]] = results(dds, contrast=c('cell_cond','THP1_PMA+LPS_RNA','THP1_PMA_RNA'), alpha=0.01)
res_diff[[3]] = results(dds, contrast=c('cell_cond','U937_PMA_RNA','U937_CTR_RNA'), alpha=0.01)
res_diff[[4]] = results(dds, contrast=c('cell_cond','U937_PMA+LPS_RNA','U937_PMA_RNA'), alpha=0.01)

# get fpkms
rna_dat_fpkm = read_tsv("inst/extdata/rna/E-MTAB-5191.genes.fpkm.htseq2.tsv")
my_samples = list(
  c("THP1 BR1 CTR_RNA","THP1 BR2 CTR_RNA","THP1 BR1 PMA_RNA","THP1 BR2 PMA_RNA"),
  c("THP1 BR1 PMA_RNA","THP1 BR2 PMA_RNA","THP1 BR1 PMA+LPS_RNA","THP1 BR2 PMA+LPS_RNA"),
  c("U937 BR1 CTR_RNA","U937 BR2 CTR_RNA","U937 BR1 PMA_RNA","U937 BR2 PMA_RNA"),
  c("U937 BR1 PMA_RNA","U937 BR2 PMA_RNA","U937 BR1 PMA+LPS_RNA","U937 BR2 PMA+LPS_RNA")
)

pos_genes = vector("list", 4)
names(pos_genes) = c("THP-1 PMA","THP-1 PMA+LPS","U937 PMA","U937 PMA+LPS")
neg_genes = pos_genes

for(i in 1:length(res_diff)) {
  
  print(i)
  
  res = tbl_df(res_diff[[i]]) # 65217 rows
  res$ensembl = rna_dat$`Gene ID`
  res$symbol = rownames_symbol
  res_filt = res %>% filter(padj<1e-8, log2FoldChange>0) %>% arrange(desc(log2FoldChange))
  
  # get genes that have <1 fpkms in control reps
  my_ids = rna_annot$`Comment[ENA_RUN]`[match(my_samples[[i]], rna_annot$`Source Name`)]
  rna_dat_fpkm_filt = rna_dat_fpkm[,match(my_ids, names(rna_dat_fpkm))] # select samples
  names(rna_dat_fpkm_filt) = my_samples[[i]]
  rna_dat_fpkm_filt = data.frame(Gene=rna_dat_fpkm$`Gene ID`, rna_dat_fpkm_filt) # add gene id
  rna_dat_fpkm_filt_baseline = rna_dat_fpkm_filt[(rna_dat_fpkm_filt[,2]<1 & rna_dat_fpkm_filt[,3]<1 & rna_dat_fpkm_filt[,4]>1 & rna_dat_fpkm_filt[,5]>1),] # filter
  
  res_filt = res_filt[res_filt$ensembl %in% rna_dat_fpkm_filt_baseline$Gene,]
  qplot(
    apply(rna_dat_fpkm_filt[match(res_filt$ensembl, rna_dat_fpkm_filt$Gene),2:3], 1, mean),
    apply(rna_dat_fpkm_filt[match(res_filt$ensembl, rna_dat_fpkm_filt$Gene),4:5], 1, mean)
  ) + xlab("Pre-Stim") + ylab("Post-Stim") + theme_thesis(25) + geom_hline(yintercept=1, color=2)
  
  # get negative set
  neg_gene = unique(rna_dat_fpkm_filt[(rna_dat_fpkm_filt[,2]<1 & rna_dat_fpkm_filt[,3]<1 & rna_dat_fpkm_filt[,4]<1 & rna_dat_fpkm_filt[,5]<1),'Gene'])
  neg_gene = mapping$hgnc_symbol[match(neg_gene, mapping$ensembl_gene_id)]
  
  # check values  
  # test_gene = rna_dat_fpkm_filt$Gene[rna_dat_fpkm_filt[,2]<1 & rna_dat_fpkm_filt[,3]<1 & rna_dat_fpkm_filt$Gene %in% res_filt$ensembl][1]
  # filter(res_filt, ensembl==test_gene)
  # filter(rna_dat_fpkm_filt, Gene==test_gene)
  
  stim_genes = unique(res_filt$symbol)
  
  # match rna response to ml matrix
  ml_ix = which(str_replace(names(my_data), "_BR[12]", "") == names(res_diff)[i])
  for(j in 1:length(ml_ix)) {
    my_data[[ml_ix[j]]]$rna = 0
    my_data[[ml_ix[j]]]$rna[my_data[[ml_ix[j]]]$Gene %in% stim_genes] = 1
    print(my_data[[ml_ix[j]]] %>% group_by(Gene) %>% summarise(Status=mean(rna)) %>% dplyr::select(Status) %>% table())
  }
  
  pos_genes[[i]] = unique(gene_list_all$hgnc_symbol[gene_list_all$hgnc_symbol %in% stim_genes])
  neg_genes[[i]] = unique(neg_gene[!neg_gene %in% stim_genes])
  
}

require(eulerr)
plot(euler(pos_genes), quantities=TRUE)


# ANALYSIS ----------------------------------------------------------------

# train on br1 and test on br2
names(my_data)
my_data_roots = unique(str_replace(names(my_data), "_BR[12]", ""))
models = vector("list", length(my_data_roots))
names(models) = my_data_roots

model_data = models
sample_ixs = models

for(i in 1:length(my_data_roots)) {
  
  # normalize?
  print(names(my_data)[which(str_replace(names(my_data), "_BR[12]", "") == my_data_roots[i] & grepl("BR1", names(my_data)))])
  dat = my_data[[which(str_replace(names(my_data), "_BR[12]", "") == my_data_roots[i] & grepl("BR1", names(my_data)))]]
  dat_norm = dat
  # dat_norm[,5:9] = normalizeQuantiles(dat[,c(5:9)])
  
  dat_trans = matrix(NA, nrow=dim(dat_norm)[1]/100, ncol=501)
  
  c_ix = 1
  for(j in 1:dim(dat_trans)[1]) {
    dat_trans[j,1:500] = unlist(dat_norm[c_ix:(c_ix+99),c(5:9)])
    dat_trans[j,501] = dat_norm$rna[c_ix]
    c_ix = c_ix+100
  }
  
  dat_trans = tbl_df(dat_trans)
  names(dat_trans)[1:500] = as.character(sapply(names(dat)[5:9], function(x) paste(x, 1:100, sep="_")))
  names(dat_trans)[501] = "Y"
  dat_trans$Y = factor(dat_trans$Y)
  
  # pheatmap(matrix(as.numeric(dat_trans[sample(which(dat_trans$Y==1),1),-dim(dat_trans)[2]]), nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))
  
  # sample from genome data
  # should the design be balanced or unbalanced?
  # https://machinelearningmastery.com/tactics-to-combat-imbalanced-classes-in-your-machine-learning-dataset/
  
  
  # CHECK FEATURE STABILITY - UNDER-SAMPLING --------------------------------
  
  fv_mat = matrix(0, nrow=500, ncol=10)
  colnames(fv_mat) = paste("Rep", 1:10)
  
  for(j in 1:10) {
    
    dat_trans_sample = dat_trans[c(which(dat_trans$Y==1), sample(which(dat_trans$Y==0 & gene_list_all$hgnc_symbol %in% neg_genes[[i]]), length(pos_genes[[i]]), replace=FALSE)),]

    task = makeClassifTask(data=dat_trans_sample, target="Y")
    lrn = makeLearner("classif.randomForest", predict.type="prob", fix.factors.prediction=TRUE)
    
    train_set = sort(c(sample(which(dat_trans_sample$Y==1),round(length(which(dat_trans_sample$Y==1))*0.8),replace=FALSE), sample(which(dat_trans_sample$Y==0),round(length(which(dat_trans_sample$Y==0))*0.8),replace=FALSE)))
    test_set = c(1:dim(dat_trans_sample)[1])[-train_set]
    
    model = train(lrn, task, subset=train_set)
    pred = predict(model, task=task, subset=test_set)
    print(performance(pred, measures=list(mmce, acc)))
    fv = generateFilterValuesData(task, method="information.gain")
    fv_mat[,j] = fv$data$information.gain
    
  }
  
  png(paste0("tmp/image_",i,"_1.png"), height=600, width=1000)
  pheatmap(fv_mat, cluster_rows=FALSE, cluster_cols=FALSE)
  dev.off()
  
  # CHECK FEATURE STABILITY - UNDER-SAMPLING --------------------------------
  
  
  # sample negative data
  sample_ixs[[i]] = c(which(dat_trans$Y==1), sample(which(dat_trans$Y==0 & gene_list_all$hgnc_symbol %in% neg_genes[[i]]), length(pos_genes[[i]]), replace=FALSE))
  dat_trans_sample = dat_trans[sample_ixs[[i]],]
  model_data[[i]] = dat_trans_sample
  
  # set up model
  task = makeClassifTask(data=dat_trans_sample, target="Y")
  lrn = makeLearner("classif.randomForest", predict.type="prob", fix.factors.prediction=TRUE)
  
  # define training/testing sets
  train_set = sort(c(sample(which(dat_trans_sample$Y==1),round(length(which(dat_trans_sample$Y==1))*0.8),replace=FALSE), sample(which(dat_trans_sample$Y==0),round(length(which(dat_trans_sample$Y==0))*0.8),replace=FALSE)))
  test_set = c(1:dim(dat_trans_sample)[1])[-train_set]
  
  # train
  model = train(lrn, task, subset=train_set)
  models[[i]] = model
  pred = predict(model, task=task, subset=test_set)
  print(performance(pred, measures=list(mmce, acc)))
  
  write_tsv(data.frame(
    gene=c(
      gene_list_all$hgnc_symbol[sample_ixs[[i]][head(pred$data[order(pred$data$prob.1, decreasing=TRUE),'id'],10)]],
      gene_list_all$hgnc_symbol[sample_ixs[[i]][tail(pred$data[order(pred$data$prob.1, decreasing=TRUE),'id'],10)]]
    )), paste0("data_",i,".txt"))
  
  # plot scores
  png(paste0("tmp/image_",i,"_2.png"))
  print(ggplot(data.frame(pred), aes(truth, prob.1)) + geom_boxplot() + theme_thesis() + xlab("Expression") + ylab("Score"))
  dev.off()
  
  # plot information gain
  fv = generateFilterValuesData(task, method="information.gain")
  png(paste0("tmp/image_",i,"_3.png"), height=600, width=1000)
  pheatmap(matrix(fv$data$information.gain, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE)
  dev.off()
  
  # av_pos = unlist(apply(filter(dat_trans_sample, Y==1) %>% dplyr::select(-Y), 2, mean))
  # av_neg = unlist(apply(filter(dat_trans_sample, Y==0) %>% dplyr::select(-Y), 2, mean))
  # av_range = range(c(av_pos, av_neg))
  # pheatmap(matrix(av_pos, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=av_range[1], to=av_range[2], by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))
  # pheatmap(matrix(av_neg, nrow=5, byrow=TRUE), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=av_range[1], to=av_range[2], by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))
  
}


# CHECK MODELS ACROSS CELL LINES ------------------------------------------

# does the thp1 model predict u937 expression?

task = makeClassifTask(data=model_data$`THP-1_Baseline`, target="Y")
pred = predict(models$U937_Baseline, task=task)
print(performance(pred, measures=list(mmce, acc)))


# HOW DOES SIGNAL DATA PERFORM? -------------------------------------------

load("data/dat_all.RData")

# take the gene sets for each model
# create an input vector with 3 different signal metrics for each data type
# check performance against binning counts

for(i in 1:4) {
  
  s_ix = which(rownames(dat_all[[1]][[1]]$res) == str_replace(my_data_roots[i], "_", "_BR1_"))
  pos_ix = match(pos_genes[[i]], gene_list_all$hgnc_symbol)
  neg_genes_filt = neg_genes[[i]][neg_genes[[i]] %in% gene_list_all$hgnc_symbol]
  neg_ix = match(sample(neg_genes_filt, length(pos_genes[[i]]), replace=FALSE), gene_list_all$hgnc_symbol)
  
  input_mat = matrix(0, nrow=length(pos_ix)+length(neg_ix), ncol=15)
  colnames(input_mat) = as.character(sapply(names(dat_all[[1]])[1:5], function(x) paste(x, names(dat_all)[1:3], sep="_")))
  rownames(input_mat) = gene_list_all$hgnc_symbol[c(pos_ix, neg_ix)]
  
  c_ix = 1
  for(j in 1:5) { # cycle through data types
    for(k in 1:3) { # cycle through collapse-to-gene methods
      input_mat[,c_ix] = dat_all[[k]][[j]]$res[s_ix,c(pos_ix, neg_ix)]
      c_ix = c_ix+1
    }
  }
  
  # input_mat = log(input_mat)
  input_mat = as.data.frame(input_mat)
  input_mat$Y = factor(c(rep(1,length(pos_ix)), rep(0,length(neg_ix))))
  
  input_mat = input_mat[!apply(input_mat, 1, function(x) any(is.na(x))),]
  
  train_set = sort(c(sample(which(input_mat$Y==1), round(length(which(input_mat$Y==1))*0.8), replace=FALSE), sample(which(input_mat$Y==0), round(length(which(input_mat$Y==0))*0.8), replace=FALSE)))
  test_set = c(1:dim(input_mat)[1])[-train_set]
  
  task = makeClassifTask(data=input_mat, target="Y")
  model = train(lrn, task, subset=train_set)
  pred = predict(model, task=task, subset=test_set)
  png(paste0("image_auc_",i,"_1.png"))
  print(ggplot(data.frame(pred), aes(truth, prob.1)) + geom_boxplot() + theme_thesis() + xlab("Expression") + ylab("Score"))
  dev.off()
  performance(pred, measures=list(mmce, acc))
  
  fv = generateFilterValuesData(task, method="information.gain")
  # fv = generateFeatureImportanceData(task, learner=lrn)
  png(paste0("image_auc_",i,"_2.png"))
  print(ggplot(fv$data, aes(name, information.gain)) + geom_bar(stat="identity") + theme_thesis(20) + xlab("") + ylab("Information"))
  dev.off()
  
}


# NEXT STEPS --------------------------------------------------------------

# try different predictors



