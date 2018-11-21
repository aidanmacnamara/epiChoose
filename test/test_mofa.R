
load("data/dat_all.RData")


# DATA PREP ---------------------------------------------------------------

# list of matrices: features as rows, samples as columns

mofa_input = vector("list", length(dat_all$max))
names(mofa_input) = names(dat_all$max)
mofa_input$H3K27ac = t(as.matrix(dat_all$max$H3K27ac$res))
mofa_input$H3K4me3 = t(as.matrix(dat_all$max$H3K4me3$res))
mofa_input$H3K27me3 = t(as.matrix(dat_all$max$H3K27me3$res))
mofa_input$ATAC = t(as.matrix(dat_all$closest$ATAC$res))
mofa_input$CTCF = t(as.matrix(dat_all$closest$CTCF$res))
mofa_input$RNA = t(as.matrix(dat_all$max$RNA$res))

# remove features with all na
for(i in 1:length(mofa_input)) {
  mofa_input[[i]] = mofa_input[[i]][apply(mofa_input[[i]], 1, function(x) !all(is.na(x))),]
}

mofa = createMOFAobject(mofa_input)
mofa

plotTilesData(mofa)

# options
# data are centred (0 mean)
# unscaled (not unit variance)
# 25 factors default
# likelihood models for each view guessed - how to check?
mofa_model@ModelOptions$likelihoods
# sparsity is used (default)

# training options
getDefaultTrainOptions()
mofa_train = prepareMOFA(mofa)
mofa_model = runMOFA(mofa_train, outfile=tempfile())
mofa_model


# VISUALISATION -----------------------------------------------------------

plotVarianceExplained(mofa_model)

plotWeightsHeatmap(
  mofa_model, 
  view = "H3K27ac", 
  factors = 1:5,
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

ggplot(filter(y, (project!="DEEP")), aes(x=x, y=y, color=project)) + geom_point(size=5, shape=17) + xlab("") + ylab("") + theme_thesis() + geom_text_repel(aes(label=annot_1), fontface="bold", size=3, force=0.5)

clusters <- clusterSamples(mofa_model, k=2, factors=1)
plotFactorScatter(mofa_model, factors=1:2, color_by=clusters)


