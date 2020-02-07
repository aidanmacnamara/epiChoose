
require(tidyverse)
require(MOFA)
require(vsn)
require(BiocParallel)
require(MultiAssayExperiment)
require(reticulate)
require(epiChoose)
require(ggrepel)

use_condaenv("general", required=TRUE)
# use_python(python="/usr/local/bin/python3", required=TRUE)
load("data/dat_all.RData")


# DATA PREP ---------------------------------------------------------------

# use row_ix and col_data from 'test/test_figure_generate.r'
# list of matrices: features as rows, samples as columns

col_data = dat_all$tss$H3K27ac$annot
col_data = data.frame(dplyr::select(col_data, Label))
row_ix = which(!is.na(col_data$Label))

inputs = list(
  t(as.matrix(dat_all$max$H3K27ac$res[row_ix,])),
  t(as.matrix(dat_all$max$H3K4me3$res[row_ix,])),
  t(as.matrix(dat_all$max$H3K27me3$res[row_ix,]))
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
col_data$rep = str_extract(col_data$Label, "BR[12]")
col_data$condition = str_extract(col_data$Label, "[[:alnum:]+]+$")
col_data$cell_type = str_extract(col_data$Label, "^[[:alnum:]-]+")
primary_ix = c(1:95,178:221)
col_data$rep[primary_ix] = "BR1"
col_data$cell_type[primary_ix] = "primary"
col_data$condition[c(1:6,37,40,41)] = "primary_monocyte"
col_data$condition[c(9,11,12,38,39)] = "primary_macrophage"
col_data$condition[c(7,8,10)] = "primary_macrophage_inflamm"
col_data$source = c(rep("Blueprint",12), rep("GSK",24),"Encode", rep("Wang",4))
col_data$group = paste(col_data$cell_type, col_data$condition, sep="_")
col_data$group[primary_ix] = col_data$condition[primary_ix]
col_data[,2:6] = lapply(col_data[,2:6], factor)
col_data = tbl_df(col_data)
col_data

annot = data.frame(col_data)
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
mofa_train@ModelOptions$numFactors = 15
mofa_train@ModelOptions
mofa_model = runMOFA(mofa_train, outfile=tempfile())
mofa_model
save(mofa_model, file="tmp/mofa_model.RData")


# VISUALISATION -----------------------------------------------------------

plotVarianceExplained(mofa_model)

plotWeightsHeatmap(
  mofa_model, 
  view = "H3K27ac", 
  factors = 1:4,
  show_colnames = FALSE
)

mofa_factors <- getFactors(mofa_model, factors=1:5, as.data.frame=FALSE)
head(mofa_factors)

# pca_res <- prcomp(mofa_factors, scale=TRUE, center=TRUE)
# pca_res_summary = summary(pca_res)
# y = data.frame(pca_res$x[,1:2])

y = data.frame(mofa_factors[,1:2])

names(y) = c("x","y")
y$annot_1 = rownames(mofa_factors)
y$project = dat_all$max$H3K27ac$annot$Project[match(y$annot_1, dat_all$max$H3K27ac$annot$Label)]

# ggplot(y, aes(x=x, y=y)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=2, force=0.5)

ggplot(y[row_ix,], aes(x=x, y=y, color=project)) + geom_point(size=5, shape=17) + xlab("") + ylab("") + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=3, force=0.5)

clusters <- clusterSamples(mofa_model, k=3, factors=1)
plotFactorScatter(mofa_model, factors=1:2, color_by=clusters)


