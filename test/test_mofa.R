
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
col_data = col_data %>% filter(!is.na(Label))

donors = list(
  monocyte = c("C0010K","C00408","C000S5","C0011I","C004SQ","C001UY"),
  macrophage = c("S0022I","S00390","S001S7"),
  inflamm = c("S0022I","S001MJ","S001S7")
)

# groups
b_ix = 1:52
c_ix = 83:163
s_ix = 53:82
w_ix = 186:189
e_ix = 164:180
n_ix = 17:24
d_ix = 181:185

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

# variance explained

plot_variance_explained = function(object, cluster=TRUE, ...) {
  
  r2_list = calculateVarianceExplained(object, ...)
  fvar_m = r2_list$R2Total
  fvar_mk = r2_list$R2PerFactor
  fvar_mk_df = reshape2::melt(fvar_mk, varnames=c("factor","view"))
  fvar_mk_df$factor = factor(fvar_mk_df$factor)
  
  if(cluster & ncol(fvar_mk) > 1) {
    hc = hclust(dist(t(fvar_mk)))
    fvar_mk_df$view = factor(fvar_mk_df$view, levels=colnames(fvar_mk)[hc$order])
  }
  
  hm = ggplot(fvar_mk_df, aes_string(x="view", y="factor")) + geom_tile(aes_string(fill="value"), color="black") + guides(fill=guide_colorbar("r2")) + scale_fill_gradientn(trans="sqrt", colors=c("gray97","darkblue"), guide="colorbar") + ylab("Latent factor") + theme(plot.title=element_text(size=17, hjust=0.5), axis.title.x=element_blank(), axis.text.x=element_text(size =11, angle=60, hjust=1, vjust=1, color="black"), axis.text.y=element_text(size=12, color="black"), axis.title.y=element_text(size=15), axis.line=element_blank(), axis.ticks=element_blank(), panel.background=element_blank())
  
  hm = hm + ggtitle("Variance explained per factor") + guides(fill=guide_colorbar("r2"))
  
  fvar_m_df = data.frame(view=factor(names(fvar_m), levels=names(fvar_m)), r2=fvar_m)
  if (cluster==TRUE & ncol(fvar_mk) > 1) {
    fvar_m_df$view = factor(fvar_m_df$view, levels=colnames(fvar_mk)[hc$order])
  }
  bplt = ggplot(fvar_m_df, aes_string(x="view", y="r2")) + ggtitle("Total variance explained per view") + geom_bar(stat="identity", fill="deepskyblue4", width=0.9) + xlab("") + ylab("r2") + scale_y_continuous(expand=c(0.01, 0.01)) + theme(plot.margin=unit(c(1,2.4,0,0), "cm"), panel.background=element_blank(), plot.title=element_text(size=17, hjust=0.5), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=12, color="black"), axis.title.y=element_text(size=13, color="black"), axis.line=element_line(size=rel(1), color="black"))
  p = plot_grid(bplt, hm, align="v", nrow=2, rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
  
}

plotFactorCor(mofa_model)
plot_variance_explained(mofa_model)
ggsave("plots/mofa_varexp.pdf",width=4,height=8)


## factor scatters
```{r,fig.width=10,fig.height=10, error=FALSE, message=FALSE, warning=FALSE}
# LF 1 and 8 separate stimulations
plotFactorScatters(model, color_by="stimulation")

# LF 2 and 7 separate celltypes 
plotFactorScatters(model, color_by="celltype")

# no donor-separatino (since donor effect was removed)
plotFactorScatters(model, color_by="donor")

# LF 6 separates timepoints
plotFactorScatters(model, color_by="timepoint")
```

# correlation network
## build correlation graph
```{r}

# useful for annotation
col_data = colData(model@InputData) %>% as_tibble(rownames="condition")

# pairwise correlations based on latent factors
corrs = getFactors(model) %>%
  t %>%
  cor(.,.,use="pairwise.complete.obs")
corrs[upper.tri(corrs,diag=TRUE)] = NA
corrs = corrs %>%
  as_tibble(rownames="cond1") %>%
  gather(cond2, cor, -cond1) %>%
  filter(!is.na(cor))

# edges bases on correlation
edges = corrs %>%
  select(cond1, cond2, cor) %>%
  mutate(distance = 1-cor) %>% 
  filter(cor>0) %>% 
  arrange(desc(cor))

# build graph from edges
graph_raw = as_tbl_graph(edges, directed=FALSE) %>%
  left_join(col_data, by=c("name"="condition"))
```

## filter graph
```{r}
graph_raw %>% 
  activate(edges) %>% 
  pull(cor) %>% 
  qplot(bins=100)

min_cor_threshold = .3

# filter correlations by established threshold(s)
graph = graph_raw %>% 
  activate(edges) %>% 
  filter(cor>min_cor_threshold) %>% 
  # gives better viz than only cor;
  # the .1 is a "pseudo-count" that no edge 
  # has weight 0 in the cluster detection
  mutate(weight=cor-min_cor_threshold+.1)
```


## plots
```{r, fig.width=20, fig.height=20}

# color = stimulation, shape = celltype
set.seed(8)
graph %>%
  ggraph(layout="fr") +
  geom_edge_link(aes(edge_width=cor),edge_alpha=.4, edge_color="gray60") +
  geom_node_point(aes(shape=celltype,color=stimulation),size=10) +
  scale_color_jco(name="stimulation") +
  theme(
    legend.justification = c(0, 0), 
    legend.position = c(0, 0)
  )
ggsave("plots/network_stimulation.pdf",width=25,height=13)

# color = donor, shape = celltype
set.seed(8)
graph %>%
  ggraph(layout="fr") +
  geom_edge_link(aes(edge_width=cor),edge_alpha=.4, edge_color="gray60") +
  geom_node_point(aes(shape=celltype,color=donor),size=10) +
  scale_color_jco(name="donor") +
  theme(
    legend.justification = c(0, 0), 
    legend.position = c(0, 0)
  ) 
ggsave("plots/network_donor.pdf",width=25,height=13)

```

# latent factors
## heatmap
```{r,fig.height=7}

mat = getFactors(model)

# make short rownames
mat = mat %>%  
  as_tibble(rownames="condition") %>% 
  left_join(
    graph %>% 
      activate(nodes) %>% 
      as_tibble(),
    by=c("condition"="name")
  ) %>% 
  select(condition,starts_with("LF")) %>%
  arrange(condition) %>% 
  column_to_rownames("condition") %>% 
  as.matrix() %>% 
  {.[,!colAlls(is.na(.))]}

paletteLength = 50
colors = colorRampPalette(c("royalblue3", "white", "red3"))(paletteLength)
breaks = c(seq(min(mat,na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
            seq(max(mat,na.rm=T)/paletteLength, max(mat,na.rm=T), length.out=floor(paletteLength/2)))

pheatmap::pheatmap(mat, color=colors, breaks=breaks)
```

## factor plots
### preparations
```{r}

# plots the values of a latent factor in the graph
plot_lf_in_graph = function(graph, lf, model) {
  
  if(!lf %in% colnames(getFactors(model))) {
    stop("Latent factor not present in the model.")
  }
  
  set.seed(1)
  graph %>% 
    activate(nodes) %>% 
    left_join(
      getFactors(model, as.data.frame=TRUE) %>% 
        filter(factor==lf),
      by=c("name"="sample")
    ) %>% 
    ggraph(layout="fr") +
    geom_edge_link(aes(edge_width=cor),edge_alpha=.4, edge_color="gray60") +
    geom_node_point(aes(shape=celltype,color=value),size=8) +
    geom_node_text(aes(label=stimulation),size=4,fontface="bold") +
    scale_color_gradient2(low="royalblue3", mid="white", high="red3", midpoint = 0, name="LF value") +
    scale_edge_width(range=c(.1,1.5)) +
    theme(
      legend.justification = c(0, 0), 
      legend.position = c(0, 0)
    ) 
}

# plot latent factor weights of a MOFA model
plot_weights = function(model,view,lf,highlight=30) {
  df = getWeights(model,view,as.data.frame=TRUE) %>% 
    filter(factor==lf) %>% 
    mutate(highlight2=ifelse(rank(-abs(value))<=highlight,TRUE,FALSE)) %>% 
    mutate(rank_value=rank(value))
  
  df %>% 
    ggplot(aes(rank_value,value,color=highlight2,label=feature)) +
    geom_point(size=1) +
    scale_color_manual(values=c("black","gray50"),limits=c(TRUE,FALSE)) +
    theme(legend.position = "none") +
    geom_text_repel(data=filter(df,highlight2))+
    ggtitle(view)
}

# beeswarm plot of factor value with annotation of samples
plot_factor_samples = function(model, lf) {
  
  df = getFactors(model, lf, as.data.frame = TRUE) %>% 
    left_join(col_data, by=c("sample"="condition")) 
  
  ggpubr::ggbarplot(df, x="sample", y="value", fill = "stimulation", color="white", 
                    sort.val = "asc", sort.by.groups = FALSE, palette="jco") +
    coord_flip() +
    theme(axis.text.y = element_text(size=8))
}

model_lf_plot = model
db = AnnotationHub::query(AnnotationHub(),c("org.Hs.eg.db"))[[1]]

# metabolite IDs to names
metabolite_anno = rowData(mae[["Metabolomics"]]) %>% 
  as_tibble(rownames="rowname") %>% 
  select(rowname, First_metabolite_name_HMDB, Ion_mz) %>% 
  mutate(name_short=ifelse(First_metabolite_name_HMDB=="",round(Ion_mz,3),First_metabolite_name_HMDB)) %>% 
  column_to_rownames("rowname")

# informative gene annotation
featureNames(model_lf_plot)$RNASeq = featureNames(model_lf_plot)$RNASeq %>% 
  {mapIds(db,.,"SYMBOL","ENSEMBL",multiVals="first")}

# informative metabolite annotation
featureNames(model_lf_plot)$Metabolomics = metabolite_anno[featureNames(model_lf_plot)$Metabolomics,]$name_short

# informative protein annotation
featureNames(model_lf_plot)$ExpressionProteomics = featureNames(model_lf_plot)$ExpressionProteomics %>% 
  {mapIds(db,.,"SYMBOL","UNIPROT",multiVals="first")}
```

## plots
```{r,fig.width=25, fig.height=15}

plot_factor = function(model, lf, graph) {
  plot_grid(
    plot_grid(
      plot_lf_in_graph(graph, lf, model),
      plot_factor_samples(model,lf), nrow=1
    ),
    plot_grid(
      plot_weights(model, "RNASeq", lf),
      plot_weights(model, "ExpressionProteomics", lf),
      plot_weights(model, "Metabolomics", lf), nrow=1
    ), ncol=1, rel_heights = c(2,1)
  ) %>% 
    plot_grid(ggdraw()+draw_label(paste0("Latent Factor ",str_replace(lf,"LF","")),fontface="bold",size=24),.,ncol=1,rel_heights = c(.1,1))
}

pdf("factors.pdf", width=25, height=15)
getFactors(model_lf_plot) %>% 
  colnames %>% 
  map(~plot_factor(model_lf_plot, .x, graph))
dev.off()

```

# enrichment for latent factors
## factors of interest
```{r, fig.width=6,fig.height=6}

# stimulation effects
getFactors(model) %>%
  as_tibble(rownames="condition") %>%
  left_join(colData(model@InputData) %>% as_tibble(rownames="condition"), by="condition") %>% 
  ggplot(aes(x=LF1,y=LF8,color=stimulation,label=condition)) +
  geom_point() +
  scale_color_jco() +
  theme(legend.position = "bottom")
ggsave("plots/stimulation_factors.pdf",width=5,height=5)

# M1 vs M0/M2 effect
lf1_samples = getFactors(model) %>%
  as_tibble(rownames="condition") %>%
  left_join(colData(model@InputData) %>% as_tibble(rownames="condition"), by="condition") %>% 
  ggplot(aes(x="",y=LF1,color=stimulation)) +
  ggbeeswarm::geom_quasirandom() +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
lf1_samples

# M2 vs M0/M1 effect
lf8_samples = getFactors(model) %>%
  as_tibble(rownames="condition") %>%
  left_join(colData(model@InputData) %>% as_tibble(rownames="condition"), by="condition") %>% 
  ggplot(aes(x="",y=LF8,color=stimulation)) +
  ggbeeswarm::geom_quasirandom() +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
lf8_samples

# celltype effect
getFactors(model) %>%
  as_tibble(rownames="condition") %>%
  left_join(colData(model@InputData) %>% as_tibble(rownames="condition"), by="condition") %>% 
  ggplot(aes(x=LF2,y=LF7,color=celltype)) +
  geom_point() +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  geom_abline(intercept=0, slope=-.25, linetype="33")
ggsave("plots/celltype_factors.pdf",width=7,height=7)

lf2_7_samples = getFactors(model) %>%
  as_tibble(rownames="condition") %>%
  mutate(LF2_minus_LF7=LF2-.25*LF7) %>% 
  left_join(colData(model@InputData) %>% as_tibble(rownames="condition"), by="condition") %>% 
  ggplot(aes(x="",y=LF2_minus_LF7,color=celltype)) +
  ggbeeswarm::geom_quasirandom() +
  scale_color_jco() +
  theme(legend.position = "bottom")+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
lf2_7_samples
```

## preparations
```{r}
# convert MSigDB gene set collections to MOFA-usable format
gmt_to_member_matrix = function(gmt_file) {
  
  gsc_list = piano::loadGSC(gmt_file)$gsc
  
  map_dfr(names(gsc_list), ~tibble(gsc=.x,gene_name=gsc_list[[.x]])) %>% 
    spread(gene_name, gene_name) %>% 
    mutate_at(-1,~ifelse(is.na(.x),0L,1L)) %>%
    mutate(gsc=str_replace_all(gsc,"^GO_|^HALLMARK_","") %>% str_sub(end=50) %>% make.names(unique=TRUE)) %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    column_to_rownames("gsc") %>% 
    as.matrix()
}

# plots feature enrichment for a factor in a view
plot_factor_enrichment_view = function(model, lf, view, feature_sets, detailed=FALSE, ...) {
  
  enrichment_res = runEnrichmentAnalysis(model, view, feature_sets, factors=lf, ...)
  
  if(detailed) {
    res = plotEnrichmentDetailed(model_enrichment, lf, feature_sets, enrichment_res) +
      theme(axis.text.y = element_text(size=rel(.8)))
  } else {
    res = plotEnrichment(model,enrichment_res,lf) +
      theme(axis.text.y = element_text(size=rel(.8)))
  }
  if(!is.null(res)) {
    res
  } else {
    ggdraw()
  }
}

rename_NA_features = function(model) {
  featureNames(model)[viewNames(model)] = featureNames(model)[viewNames(model)] %>% 
    map(function(x) {
      x[is.na(x)] = "not_available"
      make.unique(x)
    })
  model
}

# feature sets
gene_names_to_hallmark =  gmt_to_member_matrix(gsc_files["hallmark"])
gene_IDs_to_metacore = readRDS("metacore_maps/Myeloid_genes2metacore_maps.rds") %>% 
  `rownames=`(str_sub(rownames(.),end=50))
gRanges_to_metacore = readRDS("metacore_maps/Myeloid_ATACRegions2metacore_maps.rds") %>% 
  `rownames=`(str_sub(rownames(.),end=50))

# ========================================================= 
# make model with only gene names in RNA-seq and proteomics and
# =========================================================
# names in metabolomics
model_gene_names = model
db = AnnotationHub::query(AnnotationHub(),c("EnsDb","sapiens","98"))[[1]]
featureNames(model_gene_names)$RNASeq = featureNames(model_gene_names)$RNASeq %>% 
  {mapIds(db,.,"SYMBOL","GENEID",multiVals="first")}

featureNames(model_gene_names)$ExpressionProteomics = featureNames(model_gene_names)$ExpressionProteomics %>% 
  {mapIds(db,.,"SYMBOL","UNIPROTID",multiVals="first")}


# metabolite IDs to names
metabolite_anno = rowData(mae[["Metabolomics"]]) %>% 
  as_tibble(rownames="rowname") %>% 
  select(rowname, First_metabolite_name_HMDB, Ion_mz) %>% 
  mutate(name_short=ifelse(First_metabolite_name_HMDB=="",round(Ion_mz,3),First_metabolite_name_HMDB)) %>% 
  column_to_rownames("rowname")

# informative metabolite annotation
featureNames(model_gene_names)$Metabolomics = metabolite_anno[featureNames(model_gene_names)$Metabolomics,]$name_short

model_gene_names = rename_NA_features(model_gene_names)

# =========================================================
# make model with only gene IDs in RNA-seq and proteomics
# =========================================================
model_gene_IDs = model

featureNames(model_gene_IDs)$ExpressionProteomics = featureNames(model_gene_IDs)$ExpressionProteomics %>% 
  {mapIds(db,.,"GENEID","UNIPROTID",multiVals="first")}

model_gene_IDs = rename_NA_features(model_gene_IDs)
```

## plot factors of interest
```{r}
plot_factor_enrichment = function(model_gene_names, model_gene_IDs, factor, factor_values_plot) {
  
  values_and_loadings = plot_grid(
    factor_values_plot,
    plotWeights(model_gene_names, "RNASeq", factor, nfeatures = 30) + ggtitle("RNA"),
    plotWeights(model_gene_names, "ExpressionProteomics", factor, nfeatures = 20) + ggtitle("Protein"),
    plotWeights(model_gene_names, "Metabolomics", factor, nfeatures = 20) + ggtitle("Metabolomics"),
    nrow=1
  )
  enrich_hallmark = plot_grid(
    ggdraw() + draw_label("Enrichment hallmark", size=24, angle=90),
    plot_factor_enrichment_view(model_gene_names, factor, "RNASeq", gene_names_to_hallmark, transformation="none") + ggtitle("RNA"),
    plot_factor_enrichment_view(model_gene_names, factor, "ExpressionProteomics", gene_names_to_hallmark, transformation="none") + ggtitle("Protein"),
    nrow=1, rel_widths = c(.1,1,1,1)
  )
  enrich_metacore = plot_grid(
    ggdraw() + draw_label("Enrichment metacore", size=24, angle=90),
    plot_factor_enrichment_view(model_gene_IDs, factor, "RNASeq", gene_IDs_to_metacore, transformation="none") + ggtitle("RNA"),
    plot_factor_enrichment_view(model_gene_IDs, factor, "ExpressionProteomics", gene_IDs_to_metacore, transformation="none") + ggtitle("Protein"),
    plot_factor_enrichment_view(model_gene_IDs, factor, "ATACSeq", gRanges_to_metacore, transformation="none") + ggtitle("ATAC"),
    nrow=1, rel_widths = c(.1,1,1,1)
  )
  
  plot_grid(ggdraw() + draw_label(paste0("Enrichment ",factor), size=30),values_and_loadings, enrich_hallmark, enrich_metacore, ncol=1, rel_heights = c(.1,1,1,1))
}

plot_factor_enrichment(model_gene_names, model_gene_IDs, "LF1", lf1_samples)
ggsave("plots/lf1_enrich.pdf",width=30,height=20)

plot_factor_enrichment(model_gene_names, model_gene_IDs, "LF8", lf8_samples)
ggsave("plots/lf8_enrich.pdf",width=30,height=20)


## make a linear combination of latent factors
model_LF2_7_names = model_gene_names

model_LF2_7_names@Expectations$W = model_LF2_7_names@Expectations$W %>% 
  map(~cbind(.x,"LF2_minus_LF7"=.x[,"LF2"]-.25*.x[,"LF7"]))
model_LF2_7_names@Expectations$Z = cbind(model_LF2_7_names@Expectations$Z,"LF2_minus_LF7"=model_LF2_7_names@Expectations$Z[,"LF2"]-.25*model_LF2_7_names@Expectations$Z[,"LF7"])

model_LF2_7_IDs = model_gene_IDs

model_LF2_7_IDs@Expectations$W = model_LF2_7_IDs@Expectations$W %>% 
  map(~cbind(.x,"LF2_minus_LF7"=.x[,"LF2"]-.25*.x[,"LF7"]))
model_LF2_7_IDs@Expectations$Z = cbind(model_LF2_7_IDs@Expectations$Z,"LF2_minus_LF7"=model_LF2_7_IDs@Expectations$Z[,"LF2"]-.25*model_LF2_7_IDs@Expectations$Z[,"LF7"])

plot_factor_enrichment(model_LF2_7_names, model_LF2_7_IDs, "LF2_minus_LF7", lf2_7_samples)
ggsave("plots/lf2_7_enrich.pdf",width=30,height=20)

