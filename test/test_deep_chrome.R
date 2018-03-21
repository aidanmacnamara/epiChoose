
load("data/gene_list_all.RData")


# EXPORT COORDINATES ------------------------------------------------------

tss_win = gene_list_all
start(tss_win) = tss_win$transcription_start_site - 5e3
end(tss_win) = tss_win$transcription_start_site + (5e3+1)
write_tsv(data.frame(tss_win)[,1:3], "y:/links/projects/bin_regions/tss_win.bed", col_names=FALSE)


# IMPORT AND MUNGE DATA ---------------------------------------------------

# run the code in links/projects/bin_regions

# y = read_tsv("y:/links/projects/bin_regions/otar_samples/thp_1_BR1.txt", col_names=FALSE)
# thp1_df = lapply(as.list(thp1), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/", x), col_names=FALSE))
# y_2 = read_tsv("y:/sandbox/bin_regions/otar_samples/thp_1_BR1_atac.txt", col_names=FALSE)
# y = cbind(y,y_2$X5)
# names(y) = c("Seq","Start","End","Bin","H34K27ac","H3K4me3","CTCF","H3K27me3","ATAC")
# y$Bin = rep(1:100, 18086)
# y$Gene = rep(gene_list_all$hgnc_symbol, each=100)
# saved as dc_example

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

write_csv(train_set[,c(10,4:9,11)], "tmp/train.csv", col_names=FALSE)
write_csv(test_set[,c(10,4:9,11)], "tmp/test.csv", col_names=FALSE)
write_csv(validate_set[,c(10,4:9,11)], "tmp/valid.csv", col_names=FALSE)

# upload to drip and run deepchrome


# ANALYSIS ----------------------------------------------------------------

y = read_tsv("c:/Downloads/table.txt", col_names=FALSE)
y = data.matrix(y)
y[1:5,1:5]
pheatmap(y, cluster_rows=FALSE, cluster_cols = FALSE)


# MOCK EXAMPLE ------------------------------------------------------------

# add in 3 signals across datatypes for positive
hist(log(as.numeric(unlist(dc_example[,5:9]))))
pos_data = data.frame()
neg_data = data.frame()

for(i in 1:30) { # 10 positive genes * 3
  
  dat_raw_pos = sample(as.numeric(unlist(dc_example[,5:9])), size=500, replace=TRUE)
  dat_raw_neg = sample(as.numeric(unlist(dc_example[,5:9])), size=500, replace=TRUE)
  # pheatmap(t(matrix(dat_raw, nrow=100, byrow=FALSE)), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))
  
  # add signal
  dat_raw_pos[c(50:75,315:340)] = dat_raw_pos[c(50:75,315:340)]*100
  # pheatmap(t(matrix(dat_raw, nrow=100, byrow=FALSE)), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))
  pos_data = rbind(pos_data, cbind(i, 1:100, matrix(dat_raw_pos, nrow=100, byrow=FALSE), 1))
  neg_data = rbind(neg_data, cbind(i+30, 1:100, matrix(dat_raw_neg, nrow=100, byrow=FALSE), 0))
}
names(pos_data)[1] = "V1"

write_csv(rbind(pos_data[1:1000,], neg_data[1:1000,]), "tmp/train.csv", col_names=FALSE)
write_csv(rbind(pos_data[1001:2000,], neg_data[1001:2000,]), "tmp/test.csv", col_names=FALSE)
write_csv(rbind(pos_data[2001:3000,], neg_data[2001:3000,]), "tmp/valid.csv", col_names=FALSE)


# TRY KERAS ---------------------------------------------------------------

require(keras)

dat_training = read_csv("tmp/train.csv", col_names=FALSE)
dat_testing = read_csv("tmp/test.csv", col_names=FALSE)
y_training = matrix(0, nrow=20, ncol=2)
y_training[1:10,1] = 1
y_training[11:20,2] = 1
y_testing = y_training

x_training = array(NA, dim=c(20,5,100))
x_testing = x_training
c_ix = 1
for(i in 1:dim(x_training)[1]) {
  x_training[i,,] = t(dat_training[c_ix:(c_ix+99),3:7])
  x_testing[i,,] = t(dat_testing[c_ix:(c_ix+99),3:7])
  c_ix = c_ix+100
}

pheatmap(x_training[1,,], cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

pheatmap(x_training[11,,], cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(from=0, to=100, by=10), color=colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(10))

# model 1

x_training <- array_reshape(x_training, c(nrow(x_training), 500))
x_testing <- array_reshape(x_testing, c(nrow(x_testing), 500))

model <- keras_model_sequential() 
model %>% 
  layer_dense(units=256, activation='relu', input_shape=c(500)) %>% 
  layer_dropout(rate=0.4) %>% 
  layer_dense(units=128, activation='relu') %>%
  layer_dropout(rate=0.3) %>%
  layer_dense(units=2, activation='softmax')

model %>% compile(loss='categorical_crossentropy', optimizer=optimizer_rmsprop(), metrics=c('accuracy'))
history <- model %>% fit(x_training, y_training, epochs=30, batch_size=128, validation_split=0.2)
plot(history)
model %>% evaluate(x_testing, y_testing)
model %>% predict_classes(x_testing)

# model 2 

model = keras_model_sequential()
my_v = 5 # or 100?

model %>%
  layer_conv_2d(filter=my_v, kernel_size=c(3,3), padding="same", input_shape=c(5,100)) %>%  
  layer_activation("relu") %>%  
  layer_conv_2d(filter=my_v ,kernel_size=c(3,3)) %>% layer_activation("relu") %>%
  layer_max_pooling_2d(pool_size=c(2,2)) %>%  
  layer_dropout(0.25) %>%
  layer_conv_2d(filter=my_v , kernel_size=c(3,3), padding="same") %>% layer_activation("relu") %>%
  layer_conv_2d(filter=my_v, kernel_size=c(3,3) ) %>% layer_activation("relu") %>%  
  layer_max_pooling_2d(pool_size=c(2,2)) %>%  
  layer_dropout(0.25) %>%
  layer_flatten() %>%  
  layer_dense(250) %>%  
  layer_activation("relu") %>%  
  layer_dropout(0.5) %>%  
  layer_dense(2) %>%  
  layer_activation("softmax") 


opt = optimizer_adam(lr=0.0001, decay = 1e-6 )

model %>% compile(loss="categorical_crossentropy", optimizer=opt, metrics="accuracy")
summary(model)

model %>%
  fit(x_training, y_training, batch_size=32, epochs=80, validation_data = list(x_testing, y_testing), shuffle=TRUE)


