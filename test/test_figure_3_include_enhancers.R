# look at specific comparisons
# start with thp-1 pma vs. baseline vs. sanquin primary data

load("data/total_reg.RData")
load("data/dat_all.RData")
load("data/roi_reg.RData")
load("data/gene_list_all.RData")

enhancer_ix = which(roi_reg$feature_type_name=="Enhancer" & roi_reg$SYMBOL %in% gene_list_all$hgnc_symbol & roi_reg$distanceToTSS<1e3)
roi_reg_enhancer = roi_reg[enhancer_ix]


# DDS SANQUIN -------------------------------------------------------------

row_ix = which(grepl("sanquin", rownames(total_reg$H3K27ac$res), ignore.case=TRUE) & !is.na(total_reg$H3K27ac$annot$Name))

col_data = data.frame(label=rownames(total_reg$H3K27ac$res)[row_ix])
my_time = str_extract_all(col_data$label, "[0-9]+(hr[s]*|days)")
col_data$time = NA
col_data$challenge = 0
for(i in 1:length(my_time)) {
  if(grepl("LPS_T=4hrs$", col_data$label[i]) & length(my_time[[i]])>1) {
    col_data$challenge[i] = 1
    col_data$time[i] = my_time[[i]][2]
  } else if(length(my_time[[i]])==2) {
    col_data$time[i] = my_time[[i]][2]
  } else {
    col_data$time[i] = my_time[[i]][1]
  }
}
col_data$treatment = str_replace(col_data$label, "^.*?([[:alpha:]_]+)_T.*$", "\\1")
col_data$treatment = str_replace(col_data$treatment, "RPMI_", "")
col_data$donor = str_replace(col_data$label, "^.*mono_([0-9]+).*$", "\\1")
col_data$treatment[col_data$treatment=="RPMI"] = "Naive"
rownames(col_data) = rownames(dat_all$tss$H3K27ac$res)[row_ix]

# remove challenge data
c_ix = which(col_data$challenge==1)
row_ix_filt = row_ix[-c_ix]
col_data_filt = col_data[-c_ix,]

dat_add = total_reg$H3K27ac$res[row_ix_filt,enhancer_ix]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_sanquin = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~treatment+donor+time)
rld_sanquin = vst(dds_sanquin, blind=FALSE)


# DDS CELL LINES ----------------------------------------------------------

row_ix = which(grepl("thp-1|u937", rownames(total_reg$H3K27ac$res), ignore.case=TRUE))

col_data = total_reg$H3K27ac$annot[row_ix,]
col_data_filt = data.frame(dplyr::select(col_data, Label))
col_data_filt$rep = str_extract(col_data_filt$Label, "BR[12]")
col_data_filt$condition = str_extract(col_data_filt$Label, "[[:alnum:]+]+$")
col_data_filt$cell_line = str_extract(col_data_filt$Label, "^[[:alnum:]-]+")

dat_add = total_reg$H3K27ac$res[row_ix,enhancer_ix]
dat_add[is.na(dat_add)] = 0
dat_add = apply(dat_add, 2, as.integer)
dim(dat_add)
dds_cell_lines = DESeqDataSetFromMatrix(countData=t(dat_add), colData=col_data_filt, design=~rep+condition+cell_line)
rld_cell_lines = vst(dds_cell_lines, blind=FALSE)


# FC TABLE ----------------------------------------------------------------

# sanquin

design(dds_sanquin)
dds_sanquin$time_treatment = factor(paste(dds_sanquin$time, dds_sanquin$treatment, sep="_"))
design(dds_sanquin) = formula(~time_treatment+donor)
dds_sanquin = DESeq(dds_sanquin)

trts_sanquin = as.character(unique(dds_sanquin$time_treatment)[-1])
res_sanquin = data.frame(loci=colnames(dat_add))
for(i in 1:length(trts_sanquin)) {
  res = results(dds_sanquin, contrast=c("time_treatment",trts_sanquin[i],"1hr_Attached"))
  res_sanquin = cbind(res_sanquin, res$log2FoldChange, res$padj)
}
names(res_sanquin)[-1] = paste(rep(trts_sanquin, each=2), c("fc","p"), sep="_")
res_sanquin = gather(res_sanquin, "group","score", 2:dim(res_sanquin)[2])
res_sanquin$type = str_extract(res_sanquin$group, "[a-z]+$")
res_sanquin$group = str_replace(res_sanquin$group, "_[a-z]+$", "")

# gsk

design(dds_cell_lines)
dds_cell_lines$cell_condition = factor(paste(dds_cell_lines$cell_line,dds_cell_lines$condition,sep="_"))
design(dds_cell_lines) = formula(~cell_condition)
design(dds_cell_lines)

dds_cell_lines = DESeq(dds_cell_lines)
trts_thp1 = as.character(unique(dds_cell_lines$cell_condition)[2:6]) # pick out the treatment conditions
trts_u937 = as.character(unique(dds_cell_lines$cell_condition)[8:12])
res_cell_lines = data.frame(loci=colnames(dat_add))

# add contrasts across thp-1 and u937
for(i in 1:length(trts_thp1)) {
  res = results(dds_cell_lines, contrast=c("cell_condition",trts_thp1[i],"THP-1_Baseline"))
  res_cell_lines = cbind(res_cell_lines, res$log2FoldChange, res$padj)
}

for(i in 1:length(trts_u937)) {
  res = results(dds_cell_lines, contrast=c("cell_condition",trts_u937[i],"U937_Baseline"))
  res_cell_lines = cbind(res_cell_lines, res$log2FoldChange, res$padj)
}

names(res_cell_lines)[-1] = paste(rep(c(trts_thp1,trts_u937), each=2), c("fc","p"), sep="_")
res_cell_lines = gather(res_cell_lines, "group","score", 2:dim(res_cell_lines)[2])
res_cell_lines$type = str_extract(res_cell_lines$group, "[a-z]+$")
res_cell_lines$group = str_replace(res_cell_lines$group, "_[a-z]+$", "")

res_all = rbind(res_cell_lines, res_sanquin) # fc and p values
# split fc and p values to 2 data frames
res_all_p = filter(res_all, type=="p") %>% dplyr::select(-type) %>% spread(group, score)
res_all_fc = filter(res_all, type=="fc") %>% dplyr::select(-type) %>% spread(group, score)
rownames(res_all_p) = res_all_p$loci
rownames(res_all_fc) = res_all_fc$loci
res_all_p = res_all_p[,-1]
res_all_fc = res_all_fc[,-1]
all(names(res_all_p)==names(res_all_fc))


# WHAT GENE SETS TO USE? --------------------------------------------------

# my_p = gmtPathways("tmp/c7.all.v6.2.symbols.gmt") # immune set
my_p = gmtPathways("tmp/c5.bp.v6.2.symbols.gmt") # biological processes
my_p_filter = my_p[grep("monocyte|macrophage", names(my_p), ignore.case=TRUE)]
my_p_filter = c(my_p_filter, my_p[grep("brain", names(my_p), ignore.case=TRUE)])
my_p_annot = data.frame(name=names(my_p_filter), tissue=c(rep("myeloid",20),rep("brain",14)), mean_fc=NA, sig_fcs=NA)

# my_p_filter = my_p_filter[1:20]
# my_p_annot = my_p_annot[1:20,]

lapply(my_p_filter, length)
my_p_filter_tbl = as.data.frame.matrix((table(stack(my_p_filter))))
my_p_filter_tbl = cbind(rownames(my_p_filter_tbl), my_p_filter_tbl)
rownames(my_p_filter_tbl) = NULL
upset(my_p_filter_tbl, order.by="freq", nsets=20, text.scale=1.5)


# CORRELATION BETWEEN CELL LINES AND PRIMARY ------------------------------

# store correlation between thp and sanquin
trts_cell_lines = c(trts_thp1, trts_u937, "5days_BG") # add control
res_template = matrix(NA, nrow=length(trts_cell_lines), ncol=length(trts_sanquin))
colnames(res_template) = trts_sanquin
rownames(res_template) = trts_cell_lines

# store the results for each pathway
out_comp = vector("list", length(my_p_filter))
names(out_comp) = names(my_p_filter)

out_activity = matrix(NA, nrow=length(my_p_filter), ncol=length(c(trts_cell_lines, trts_sanquin)))
colnames(out_activity) = c(trts_cell_lines, trts_sanquin)
rownames(out_activity) = names(my_p_filter)

for(i in 1:length(out_comp)) {
  
  res_copy = res_template
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]]) # get the gene indices for the pathway
  if(is_empty(g_ix)) next
  
  pathway_fc = res_all_fc[g_ix,] # get the fold changes
  pathway_p = res_all_p[g_ix,] # get the p-values
  my_p_annot$mean_fc[i] = mean(as.numeric(unlist(pathway_fc)), na.rm=TRUE) # what is the mean fc for the pathway?
  my_p_annot$sig_fcs[i] = sum(pathway_p<0.05, na.rm=TRUE) / length(unlist(pathway_p)) # how many sig fcs for the pathway?
  
  pathway_fc[pathway_p > 0.05] = NA # set non-significant (padj) fcs to na
  
  fcs_per_condition = apply(pathway_fc, 2, function(x) sum(!is.na(x))) / dim(pathway_fc)[1]
  out_activity[i,] = fcs_per_condition[match(colnames(out_activity),names(fcs_per_condition))]
  
  for(j in 1:length(trts_cell_lines)) {
    for(k in 1:length(trts_sanquin)) {
      
      x = pathway_fc[,which(names(pathway_fc)==trts_cell_lines[j])]
      y = pathway_fc[,which(names(pathway_fc)==trts_sanquin[k])]
      # plot(x,y)
      
      fisher_mat = matrix(
        c(
          sum(x>0 & y>0, na.rm=TRUE),
          sum(x<0 & y>0, na.rm=TRUE),
          sum(x>0 & y<0, na.rm=TRUE),
          sum(x<0 & y<0, na.rm=TRUE)
        ),
        nrow=2,
        dimnames=list(cell_line=c("up","down"),primary_cell=c("up","down"))
      )
      
      fisher_res = fisher.test(fisher_mat, alternative="greater")
      res_copy[j,k] = fisher_res$p.value
      
    }
  }
  
  out_comp[[i]] = res_copy

}

ggplot(my_p_annot, aes(x=tissue, y=mean_fc)) + geom_boxplot() + theme_thesis(20) + xlab("Tissue") + ylab("Mean FC per Pathway")
ggplot(my_p_annot, aes(x=tissue, y=sig_fcs)) + geom_boxplot() + theme_thesis(15) + xlab("Tissue") + ylab("Proportion significant FCs per Pathway")

# make out_comp long to look at correlation sets
out_comp_long = data.frame()
my_labels = as.character(outer(trts_cell_lines, trts_sanquin, paste, sep="_VS_"))
for(i in 1:length(out_comp)) {
  print(i)
  out_comp_long = rbind(
    out_comp_long,
    data.frame(
      `P-value` = as.numeric(out_comp[[i]]),
      Label = my_labels,
      Pathway = names(out_comp)[i], check.names=FALSE
    )
  )
}
out_comp_long$tissue = my_p_annot$tissue[match(out_comp_long$Pathway, my_p_annot$name)]
out_comp_long = tbl_df(out_comp_long)


# CONDITIONS VS PATHWAYS --------------------------------------------------

pheatmap(out_activity)


# WHAT ARE THE HITS? ------------------------------------------------------

# the overall distribution of fisher exact test p-values
out_comp_long %>% ggplot(aes(x=`P-value`)) + geom_histogram(binwidth=0.01, color="black") + theme_thesis()

# are related primary cells more highly correlated across myeloid pathways compared to brain pathways?
out_comp_long %>% filter(Label=="5days_BG_VS_6days_Naive") %>% ggplot(aes(x=tissue, y=`P-value`)) + geom_boxplot() + theme_thesis()

cutoff = 0.05
out_comp_long %>% filter(`P-value` < cutoff) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N))
out_comp_long %>% filter(`P-value` < cutoff) %>% group_by(Label) %>% summarise(N=n()) %>% arrange(desc(N))
out_comp_long %>% filter(`P-value` < cutoff, grepl("U937", Label), tissue=="myeloid")

select_ps = out_comp_long %>% filter(`P-value` < cutoff, grepl("THP-1|U937",Label)) %>% group_by(Pathway) %>% summarise(N=n()) %>% arrange(desc(N)) %>% select(Pathway) %>% unlist %>% as.character()

for(j in 1:length(select_ps)) {
  
  # i = which(names(my_p_filter)==select_ps[j])
  i = which(names(my_p_filter)=="GO_MONOCYTE_CHEMOTAXIS")
  out_comp[[i]]
  g_ix = which(rownames(res_all_fc) %in% my_p_filter[[i]])
  dat_fc = res_all_fc[g_ix,]
  pheatmap(dat_fc, main=names(my_p_filter)[i], fontsize_row=5)
  
}

gsea_in = res_all_fc$`5days_LPS`
names(gsea_in) = rownames(res_all_fc)
res_gsea <- fgsea(pathways=my_p, stats=gsea_in, nperm=1000)
res_tidy = res_gsea %>% as_tibble() %>% arrange(desc(NES))
res_tidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% DT::datatable()


# FUNCTIONS ---------------------------------------------------------------

pca_and_plot <- function(rld, annot_1, annot_2) {
  
  if(class(rld)=="matrix") {
    y = rld
  } else {
    y = t(assays(rld)[[1]])
  }
  
  dim(y)
  y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(y)
  
  pca_res <- prcomp(y, scale=TRUE, center=TRUE)
  pca_res_summary = summary(pca_res)
  yy = data.frame(pca_res$x[,1:2])
  names(yy) = c("x","y")
  yy$annot_1 = annot_1
  yy$annot_2 = annot_2
  
  if(is.na(annot_2)) {
    my_plot = ggplot(yy, aes(x=x, y=y)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")
  } else {
    my_plot = ggplot(yy, aes(x=x, y=y, color=annot_2)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=5, force=0.5) + theme(legend.position="none")
  }
  return(my_plot)
  
}


