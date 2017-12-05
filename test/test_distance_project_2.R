# start off with mask_data

load("data/dat.RData")
load("data/column_annotation/gene_list_all.RData")
marks = names(dat) # what marks to select


# GET RNA DATA ------------------------------------------------------------

u937 = list.files("z:/links/bix-analysis-stv/2016/CTTV/U937/data/Outputs/star/bam_files/genecounts/")
thp1 = list.files("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/")

u937_df = lapply(as.list(u937), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/U937/data/Outputs/star/bam_files/genecounts/", x), col_names=FALSE))
thp1_df = lapply(as.list(thp1), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/", x), col_names=FALSE))

# merge into sample per column
my_names = u937_df[[1]]$X1
r_df = data.frame(do.call("cbind", lapply(u937_df, function(x) x$X2)))
r_df = cbind(NA, my_names, r_df)
r_df = cbind(r_df, data.frame(do.call("cbind", lapply(thp1_df, function(x) x$X2[match(my_names, x$X1)]))))

# get names
sample_info = read_csv("z:/links/bix-analysis-stv/2016/CTTV/U937/data/sampleInfo.csv", col_names=FALSE)

all_names = sample_info$X1[c(match(str_extract(u937, "^[[:alnum:]]+"), sample_info$X32), match(str_extract(thp1, "^[[:alnum:]]+"), sample_info$X32))]

names(r_df) = c("NULL", "Gene Name", all_names)

# change names so they match with epigenatic labels
names(r_df) = str_replace(names(r_df), "THP1", "THP-1")
names(r_df) = str_replace(names(r_df), "CTR", "Baseline")
names(r_df) = str_replace(names(r_df), "_RNA", "")
names(r_df) = str_replace(names(r_df), "Baseline\\+LPS", "LPS")


# SUMMARISE ENSEMBL REG ---------------------------------------------------

start_data = dat

# sample labels
single_labels = rownames(start_data[[1]]$res)
# sample_ix = c(grep("monocyte", single_labels, ignore.case=TRUE), grep("macrophage", single_labels, ignore.case=TRUE), 196:219) # project 2 data + relevant monocyte/macrophage data
sample_ix = grep("THP-1.*PMA", single_labels)
single_labels = single_labels[sample_ix]
# group_labels = c(rep("Blueprint",33), rep("ENCODE",24))
# group_labels[12] = "ENCODE"

# slice matrices
for(i in 1:length(start_data)) {
  start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
  start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
}

rna_add = prep_rna(r_df, gene_list_all$hgnc_symbol, single_labels, names(r_df)[3:26])

start_data[[6]] = rna_add
names(start_data)[6] = "RNA"


# ALL DATA SUMMARY --------------------------------------------------------

df_1 = tbl_df(data.frame(Gene=colnames(start_data[[1]]$res), do.call("cbind", lapply(start_data[c(1:3,6)], function(x) t(x$res)))))
df_rep_1 = dplyr::select(df_1, Gene, contains("BR1")) %>% mutate(Rep=1)
df_rep_2 = dplyr::select(df_1, Gene, contains("BR2")) %>% mutate(Rep=2)

names(df_rep_1) = str_replace_all(names(df_rep_1), "_BR[12]", "")
names(df_rep_2) = str_replace_all(names(df_rep_2), "_BR[12]", "")

df_1 = tbl_df(rbind(df_rep_1, df_rep_2))

names(df_1) = str_replace(names(df_1), "\\.1$", "_H3K4me3")
names(df_1) = str_replace(names(df_1), "\\.2$", "_H3K27me3")
names(df_1) = str_replace(names(df_1), "\\.3$", "_RNA")

df_1[,-c(1, dim(df_1)[2])] = log(df_1[,-c(1, dim(df_1)[2])]+1)


# LPS GENE SET ------------------------------------------------------------

load("../epiView/data/msig_go_bp.RData")
m_up = read_tsv("c:/Downloads/tmp/otar_020_tmp/GSE3982_CTRL_VS_LPS_4H_MAC_UP.txt", skip=2, col_names=FALSE)
m_down = read_tsv("c:/Downloads/tmp/otar_020_tmp/GSE3982_CTRL_VS_LPS_4H_MAC_DN.txt", skip=2, col_names=FALSE)

gene_ix = which(df_1$Gene %in% unique(c(m_up$X1, m_down$X1, unlist(msig_go_bp[grep("LIPOPOLYSACCHARIDE", names(msig_go_bp))]))))
# gene_ix = which(df_1$Gene %in% unlist(m_down$X1))

df_1_filtered = df_1[gene_ix,]

ggplot(df_1_filtered, aes(x=THP.1_PMA_RNA, y=THP.1_PMA.LPS_RNA, color=factor(Rep))) + geom_point() + theme_thesis(20) + geom_abline(slope=1) # + geom_text_repel(aes(label=Gene), fontface="bold", force=0.5)

# is there any difference in pma signature between up/non-regulated genes?

df_1_filtered$rna_diff = ifelse(df_1_filtered$THP.1_PMA.LPS > df_1_filtered$THP.1_PMA, "Yes", "No")
ggplot(df_1_filtered, aes(x=rna_diff, y=THP.1_PMA)) + geom_boxplot() + theme_thesis(20)


# COMBINATIONS ------------------------------------------------------------

ggplot(df_1_filtered, aes(x=THP.1_PMA_H3K27me3/THP.1_PMA, y=THP.1_PMA.LPS_RNA/THP.1_PMA_RNA), color=Rep) + geom_point(na.rm =TRUE) + theme_thesis(20) + xlab("PMA Epigenetics") + ylab("Expression Change") # + geom_text_repel(aes(label=label), fontface="bold", force=0.5) + scale_x_continuous(limits=c(0.5,1.25))

ggplot(df_1_filtered, aes(x=THP.1_PMA_H3K4me3/THP.1_PMA_H3K27me3, y=THP.1_PMA.LPS_RNA/THP.1_PMA_RNA), color=Rep) + geom_point(na.rm =TRUE) + theme_thesis(20) + xlab("PMA Epigenetics") + ylab("Expression Change") # + geom_text_repel(aes(label=label), fontface="bold", force=0.5) + scale_x_continuous(limits=c(0.5,1.25))

ggplot(df_1_filtered, aes(x=THP.1_PMA_H3K27me3, y=THP.1_PMA.LPS_RNA/THP.1_PMA_RNA), color=Rep) + geom_point(na.rm =TRUE) + theme_thesis(20) + xlab("PMA Epigenetics") + ylab("Expression Change") # + geom_text_repel(aes(label=label), fontface="bold", force=0.5) + scale_x_continuous(limits=c(0.5,1.25))




df_1_filtered$label = df_1_filtered$Gene
df_1_filtered$label[df_1_filtered$THP.1_BR1_PMA_H3K27me3/df_1_filtered$THP.1_BR1_PMA < 0.85] = ""

ggplot(df_1_filtered, aes(x=THP.1_BR1_PMA_H3K27me3/THP.1_BR1_PMA, y=THP.1_BR1_PMA.LPS_RNA/THP.1_BR1_PMA_RNA)) + geom_point() + theme_thesis(20) + xlab("PMA Epigenetics") + ylab("Expression Change") + geom_text_repel(aes(label=label), fontface="bold", force=0.5) + scale_x_continuous(limits=c(0.5,1.25))

# are these genes primed compared to background?

mean_auc = mean(df_1_filtered$THP.1_BR1_PMA, na.rm=TRUE)
mean_auc_all = rep(NA,10000)
for(i in 1:length(mean_auc_all)) {
  s_ix = sample(1:dim(df_1)[1], length(gene_ix), replace=FALSE)
  mean_auc_all[i] = mean(df_1$THP.1_BR1_PMA[s_ix], na.rm=TRUE)
}

ggplot(data.frame(Mean=mean_auc_all), aes(Mean)) + geom_density() + theme_thesis() + geom_point(x=mean_auc, y=0, size=5, shape=17)





