# start off with merged, normalised data - chip + rna
# taken from test_masking.r

# start_data = max reg region over gene body


# LOOP THROUGH RELEVANT LUNG GO TERMS -------------------------------------

require(qusage)
msig_go_bp = read.gmt("c5.all.v6.1.symbols.gmt")

lung_go = data.frame(
  ids = c(
    "GO:0006855","GO:0042908","GO:0006281","GO:0004033","GO:0006805","GO:0009494","GO:0017144",
    "GO:0002314","GO:0002377","GO:0002467","GO:0016064","GO:0019882","GO:0042571","GO:0006915",
    "GO:0008219","GO:0000302","GO:0072593","GO:0006802","GO:0006979","GO:0034599","GO:0010193",
    "GO:0070994","GO:0000303","GO:0006801","GO:0006749","GO:0033355","GO:0034635","GO:0036347",
    "GO:1901370","GO:0006809","GO:0071731","GO:0007263","GO:0035713","GO:0046210","GO:0042744",
    "GO:0050665","GO:0045071","GO:0051607","GO:0050691","GO:0045087","GO:0045088","GO:0045089",
    "GO:0002218","GO:0002227","GO:0008063","GO:0034121","GO:0002224","GO:0001816","GO:0001817",
    "GO:0001819","GO:0002775","GO:0002794","GO:0002777","GO:0032602","GO:0032642","GO:0001875",
    "GO:0071222","GO:0038187","GO:1904019"
  ),
  names=NA,
  genes=NA
)

# add names
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

get_names = getBM(attributes=c("name_1006","go_id"), filters="go", values=lung_go, mart=ensembl)
lung_go$names = get_names$name_1006[match(lung_go$ids, get_names$go_id)]

for(i in 1:length(lung_go$ids)) {
  print(paste("Query", i))
  go_res = getBM(attributes=c("hgnc_symbol","go_id"), filters="go", values=lung_go$ids[i], mart=ensembl)
  if(!dim(go_res)[1]) next
  lung_go$genes[i] = list(unique(go_res$hgnc_symbol))
}

lung_go = tbl_df(lung_go)


# EXAMPLE -----------------------------------------------------------------

load("z:/sandbox/epiChoose/r_data/column_annotation/gene_list_all.RData")

# pick out the lung and blood samples
sample_ix = c(1:17,42:44,72)+178

# sample labels
single_labels = rownames(start_data[[1]]$res)[sample_ix]
group_labels = c(rep("GSK",19), rep("ENCODE",2))

diff_test <- function(dat, parametric=TRUE, log_s=FALSE) {
  
  if(log_s) {
    dat$Score = log(dat$Score+1)
  }
  
  if(parametric) {
    y = pairlist(pairwise.t.test(dat$Score, dat$`Cell Line`, p.adjust.method="BH", paired=TRUE))
    y_df = tidy(y[[1]])
  } else {
    y = conover.test(dat$Score, dat$`Cell Line`, method="bh", table=TRUE)
    y_df = data.frame(comps=y$comparisons, P=y$P.adjusted, test_stat=y$T)
  }
  return(y_df)
}

all_res = data.frame()
all_res_means = data.frame()

go_list = lung_go$genes
# go_list = msig_go_bp

go_names = lung_go$names
# go_names = names(msig_go_bp)

for(i in 1:length(go_list)) {
  
  if(i %in% c(5598)) next
  
  genes = go_list[[i]]
  if(is.null(genes)) next
  
  col_ix = which(gene_list_all$hgnc_symbol %in% genes)
  if(length(col_ix)<2) next
  
  end_data = start_data
  
  # slice matrices
  for(j in 1:length(end_data)) { # each data type
    end_data[[j]]$res = end_data[[j]]$res[sample_ix,col_ix]
  }
  
  res = dist_mat(end_data, comp_ix=list(c(1:14,18:21), 15), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), my_title=go_names[i], plot_res=TRUE, use_corr=TRUE, font_size=20)
  
  comp_ix = c(1,3,15)
  y_all = data.frame()
  
  for(j in 1:length(end_data)) {
    x = end_data[[j]]$res[comp_ix,]
    if(all(is.na(x))) next
    y = melt(as.matrix(x))
    y = cbind(y, names(end_data[j]))
    names(y) = c("Cell Line", "Gene", "Score", "Assay")
    if(length(which((filter(y, !is.na(Score)) %>% dplyr::select(Gene) %>% table()) == length(unique(y$`Cell Line`)))) < 2) next # if there are < 2 genes with observations from all cell lines, skip ...
    y_all = rbind(y_all, y)
  }
  y_all = tbl_df(y_all)
  # y_all = filter(y_all, !is.na(Score))
  
  ggplot(y_all, aes(x=Gene, y=Score)) + geom_bar(aes(fill=`Cell Line`), position="dodge", stat="identity") + theme_thesis(10) + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Assay, nrow=2, scales="free")
  
  ggplot(y_all, aes(x=`Cell Line`, y=log(Score+1))) + geom_boxplot() + theme_thesis(10) + facet_wrap(~Assay, nrow=2, scales="free") + theme(axis.text.x=element_text(angle=45, hjust=1))
  
  # ggplot(y_all, aes(x=`Cell Line`, y=Score)) + geom_point() + geom_line(aes(group=Gene)) + theme_thesis(10) + facet_wrap(~Assay, nrow=2, scales="free")
  
  all_res_means = rbind(all_res_means, data.frame(y_all %>% group_by(Assay, `Cell Line`) %>% summarise(mean=mean(Score, na.rm=TRUE)), go=go_names[i]))
  all_res = rbind(all_res, data.frame(y_all %>% group_by(Assay) %>% do(diff_test(.,log_s=TRUE)), go=go_names[i], n=length(go_list[[i]]))) 
  
}

all_res = tbl_df(all_res)


# ADD LFC -----------------------------------------------------------------

calc_fc <- function(x, all_res_means) {
  g1_score = all_res_means$mean[all_res_means$Assay==x[1] & all_res_means$go==x[5] & all_res_means$Cell.Line==x[2]]
  g2_score = all_res_means$mean[all_res_means$Assay==x[1] & all_res_means$go==x[5] & all_res_means$Cell.Line==x[3]]
  return(log(g1_score, base=2)-log(g2_score, base=2))
}

# add log fold change
all_res$LFC = apply(all_res, 1, calc_fc, all_res_means)


# ANALYSE -----------------------------------------------------------------

# check significant rna + h3k27ac results with a bonferroni cut-off
View(filter(all_res, Assay=="RNA" | Assay=="H3K27ac", p.value<0.05) %>% arrange(p.value))


# DISTANCE BY REGULATORY TYPE ---------------------------------------------

# look at overall pc plots splitting by segmentation categories
