
require(tidyverse)
require(MOFA)
require(MOFAdata)
require(vsn)
require(BiocParallel)
require(MultiAssayExperiment)
require(reticulate)
require(epiChoose)
require(ggrepel)
require(cowplot)
require(readxl)
require(pheatmap)

use_condaenv("base", required=TRUE)
# use_python(python="/usr/local/bin/python3", required=TRUE)
load("data/dat_all.RData")


# DATA PREP ---------------------------------------------------------------

# use row_ix and col_data from 'test/test_figure_generate.r'
# list of matrices: features as rows, samples as columns

col_data = dat_all$tss$H3K27ac$annot
col_data = data.frame(dplyr::select(col_data, Label))
col_data = col_data %>% filter(!is.na(Label))

donors = list(
  monocyte = c("C0010K","C00408","C000S5","C0011I","C004SQ","C001UY"),
  macrophage = c("S0022I","S00390","S001S7"),
  inflamm = c("S0022I","S001MJ","S001S7")
)

# groups
b_ix = 1:52 # blueprint
c_ix = 83:163 # cell lines
s_ix = 53:82 # novakovic
w_ix = 186:189 # wang
e_ix = 164:180 # encode
n_ix = 17:24 # saeed
d_ix = 181:185 # deep
t_ix = 136:147 # thp-1

# rep
col_data$rep = str_extract(col_data$Label, "BR[12]")
col_data$rep[is.na(col_data$rep)] = "BR1"

# donor
col_data$donor = NA
col_data$donor[b_ix] = str_extract(col_data$Label[b_ix], "^[[:alnum:]]+")
col_data$donor[s_ix] = str_replace(col_data$Label[s_ix], "^.*_mono_([0-9]+)_.*$", "\\1")
col_data$donor[w_ix] = str_extract(col_data$Label[w_ix], "[0-9]+")
col_data$donor[c(c_ix,e_ix)] = "None"

# condition
col_data$condition = NA
col_data$condition[c_ix] = str_extract(col_data$Label[c_ix], "[[:alnum:]+]+$")
col_data$condition[c(b_ix,e_ix)] = "Baseline"
col_data$condition[s_ix] = str_replace(col_data$Label[s_ix], "^.*?([[:alpha:]_]+)_T.*$", "\\1")
col_data$condition = str_replace(col_data$condition, "RPMI_", "")
col_data$condition[col_data$condition=="RPMI"] = "Baseline"
col_data$condition[col_data$condition=="Attached"] = "Baseline"
col_data$condition[n_ix] = c("BG", "LPS", "Baseline", "Baseline", "BG", "LPS", "Baseline", "Baseline")
col_data$condition[w_ix] = "Baseline"

# cell type
col_data$cell_type = NA
col_data$cell_type[c(c_ix,e_ix)] = str_extract(col_data$Label[c(c_ix,e_ix)], "^[[:alnum:]-]+")
col_data$cell_type[b_ix] = str_extract(col_data$Label[b_ix], "[[:alnum:]\\s-\\(\\)]+$")
col_data$cell_type[b_ix] = str_replace(col_data$cell_type[b_ix], "^\\s+", "")
col_data$cell_type[c(s_ix,n_ix)] = "monocyte"
col_data$cell_type[c(s_ix,n_ix)][col_data$time[c(s_ix,n_ix)]!="0days"] = "macrophage"
col_data$cell_type[w_ix] = str_extract(col_data$Label[w_ix], "Monocyte|Macrophage")
col_data$cell_type[w_ix] = tolower(col_data$cell_type[w_ix])

# time
my_time = str_extract_all(col_data$Label, "[0-9]+(hr[s]*|days)")
col_data$time = NA
for(i in 1:length(my_time)) {
  if(grepl("LPS_T=4hrs$", col_data$Label[i]) & length(my_time[[i]])>1) {
    col_data$time[i] = my_time[[i]][2]
  } else if(length(my_time[[i]])==2) {
    col_data$time[i] = my_time[[i]][2]
  } else {
    col_data$time[i] = my_time[[i]][1]
  }
}

# source
col_data$source = NA
col_data$source[b_ix] = "blueprint"
col_data$source[c_ix] = "gsk"
col_data$source[n_ix] = "saeed"
col_data$source[s_ix] = "novakovic"
col_data$source[w_ix] = "wang"
col_data$source[e_ix] = "encode"

# group
col_data$group = paste(col_data$cell_type, col_data$condition, sep="_")

col_data[,2:8] = lapply(col_data[,2:8], factor)
col_data = tbl_df(col_data)
col_data
col_data_filt = col_data[-d_ix,]

# manuscript data
dat_ix = which(
  grepl("thp-1|u937", col_data_filt$Label, ignore.case=TRUE) | (col_data_filt$donor %in% unlist(donors))
)
dat_ix = c(dat_ix, which(col_data_filt$Label=="Monocytes-CD14+_Broad"))
dat_ix = c(dat_ix, which(col_data_filt$source=="wang"))
col_data_filt[dat_ix,]

inputs = list(
  t(as.matrix(dat_all$max$H3K27ac$res[match(col_data_filt$Label, rownames(dat_all$max$H3K27ac$res)),])),
  t(as.matrix(dat_all$max$H3K4me3$res[match(col_data_filt$Label, rownames(dat_all$max$H3K27ac$res)),])),
  t(as.matrix(dat_all$max$H3K27me3$res[match(col_data_filt$Label, rownames(dat_all$max$H3K27ac$res)),]))
  # t(as.matrix(dat_all$closest$ATAC$res)),
  # t(as.matrix(dat_all$closest$CTCF$res))
  # t(as.matrix(dat_all$max$RNA$res))
)
names(inputs) = names(dat_all$max)[1:3]

for(i in 1:length(inputs)) {
  rownames(inputs[[i]]) = paste(names(inputs)[i], rownames(inputs[[i]]), sep="_")
}

lapply(inputs, dim)
lapply(inputs, function(x) x[1:5,])
lapply(inputs, limma::plotMA)

# finish column data
annot = data.frame(col_data_filt)
rownames(annot) = colnames(inputs[[1]])
mae = MultiAssayExperiment(experiments=inputs, colData=annot)

# which features have no information across all assays
features_na = sapply(1:18246, function(x) all(is.na(unlist(assays(mae[x,,])))))
table(features_na)
mae_filt = mae[!features_na,,]
mae_filt

mofa = createMOFAobject(mae_filt)
mofa
# plotDataOverview(mofa)

# training options
getDefaultTrainOptions()
mofa_train = prepareMOFA(mofa)
mofa_train@ModelOptions$numFactors = 15
mofa_train@ModelOptions
mofa_model = runMOFA(mofa_train, outfile=tempfile())
mofa_model
save(mofa_model, file="tmp/mofa_model.RData")


# MOFA SLICE --------------------------------------------------------------

# just look at thp-1 here

inputs = list(
  t(as.matrix(dat_all$max$H3K27ac$res[match(col_data$Label[t_ix], rownames(dat_all$max$H3K27ac$res)),])),
  t(as.matrix(dat_all$max$H3K4me3$res[match(col_data$Label[t_ix], rownames(dat_all$max$H3K27ac$res)),])),
  t(as.matrix(dat_all$max$H3K27me3$res[match(col_data$Label[t_ix], rownames(dat_all$max$H3K27ac$res)),]))
  # t(as.matrix(dat_all$closest$ATAC$res)),
  # t(as.matrix(dat_all$closest$CTCF$res))
  # t(as.matrix(dat_all$max$RNA$res))
)
names(inputs) = names(dat_all$max)[1:3]

for(i in 1:length(inputs)) {
  rownames(inputs[[i]]) = paste(names(inputs)[i], rownames(inputs[[i]]), sep="_")
}

lapply(inputs, dim)
lapply(inputs, function(x) x[1:5,])
lapply(inputs, limma::plotMA)

# finish column data
annot = data.frame(col_data_filt[t_ix,])
rownames(annot) = colnames(inputs[[1]])
mae = MultiAssayExperiment(experiments=inputs, colData=annot)

# which features have no information across all assays
features_na = sapply(1:18246, function(x) all(is.na(unlist(assays(mae[x,,])))))
table(features_na)
mae_filt = mae[!features_na,,]
mae_filt

mofa = createMOFAobject(mae_filt)
mofa
plotDataOverview(mofa)

# training options
getDefaultTrainOptions()
mofa_train = prepareMOFA(mofa)
mofa_train@ModelOptions$numFactors = 5
mofa_train@ModelOptions
mofa_model = runMOFA(mofa_train, outfile=tempfile())
mofa_model
save(mofa_model, file="tmp/mofa_model_slice.RData")


# VISUALISATION -----------------------------------------------------------

plotVarianceExplained(mofa_model)

plotWeightsHeatmap(
  mofa_model, 
  view = "H3K27ac", 
  factors = 1:4,
  show_colnames = FALSE
)

mofa_factors = getFactors(mofa_model, factors=1:5, as.data.frame=FALSE)
head(mofa_factors)

# pca_res = prcomp(mofa_factors, scale=TRUE, center=TRUE)
# pca_res_summary = summary(pca_res)
# y = data.frame(pca_res$x[,1:2])

y = data.frame(mofa_factors[,1:2])

names(y) = c("x","y")
y$annot_1 = col_data_filt$group
y$project = col_data_filt$source

ggplot(y[dat_ix,], aes(x=x, y=y, color=project)) + geom_point(size=5, shape=17) + xlab("") + ylab("") + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=3, force=0.5)

ggplot(y, aes(x=x, y=y, color=project)) + geom_point(size=5, shape=17) + xlab("") + ylab("") + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=3, force=0.5)

clusters = clusterSamples(mofa_model, k=3, factors=1)
plotFactorScatter(mofa_model, factors=1:2, color_by=clusters)


# MOFA EXPLORATION --------------------------------------------------------

# binary matrix with feature sets in rows and features in columns
data("reactomeGS")

# make one for go bp
bp = gmtPathways("tmp/gene_sets/c5.bp.v6.2.symbols.gmt")
bp_table = as.matrix(as.data.frame.matrix(t((table(stack(bp))))))

# mapping file
mapping = read_excel("../network_analysis/gwas/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx")

for(i in 1:length(viewNames(mofa_model))) {
  
  my_view = viewNames(mofa_model)[i]
  print(my_view)
  
  plot_data_heatmap(mofa_model=mofa_model, plot_view=my_view, my_factor=1, features=50, include_weights=FALSE, sample_ix=dat_ix)

  # map the feature names
  fs = bp_table
  colnames(fs) = paste(my_view, colnames(fs), sep="_")
  fs[1:5,1:5]
  
  # perform enrichment analysis
  gsea = runEnrichmentAnalysis(mofa_model, view=my_view, feature.sets=fs, alpha=0.01)
  
  # plot enrichment per factor
  plotEnrichmentBars(gsea, alpha=0.01)
  
  interesting_factors = 1:2
  fsea_plots <- lapply(interesting_factors, function(factor) {
    plotEnrichment(
      mofa_model,
      gsea,
      factor = factor,
      alpha = 0.01,
      max.pathways = 10 # The top number of pathways to display
    )
  })
  print(cowplot::plot_grid(fsea_plots[[1]], fsea_plots[[2]], ncol=1, labels=paste("Factor", interesting_factors)))
  
  
  
}
