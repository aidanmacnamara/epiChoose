# start off with mask_data

load("data/roi_reg.RData")
load("data/total_data.RData")
load("data/gene_list_all.RData")

# try only promoter-annotated regions
dat_max_p = total_data
for(i in 1:length(dat_max_p[1:5])) {
  print(paste("Processing data type", names(dat_max_p)[i]))
  dat_max_p[[i]]$res = convert_reg_matrix(dat_max_p[[i]]$res[,roi_reg$feature_type_name=="Promoter"], roi_reg[roi_reg$feature_type_name=="Promoter"], gene_list_all, reg_window=0, summ_method="max")
}

single_labels = rownames(total_data$H3K27ac$res)
group_labels = factor(total_data$H3K27ac$annot$Project)

plot_pca(dat_max_p$H3K27ac$res, annot_1=single_labels, annot_2=group_labels, out_file="out.png")

# look at counts
# counts of all data types across promoters?
# plan - look at counts first and then try and repeat with auc

data_gsk = read_excel("inst/extdata/data_gsk.xlsx")



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


