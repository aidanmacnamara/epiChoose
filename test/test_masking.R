
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
load_all()


# ALL DATA ----------------------------------------------------------------

load("r_data/all_data.RData")
names(all_data) = c("RNA","H3k27ac","H3K4me3","H3K27me3")

# ROI ---------------------------------------------------------------------

# the regultory regions (2kb window around tsss)
load("r_data/roi.RData")
load("r_data/gene_list_go.RData")


# ANALYSIS ----------------------------------------------------------------

group_labels = c(rep("GSK",16), rep("BLUEPRINT",34))
single_labels = rownames(all_data[[2]]$res)

# total plot
pca_data = prep_for_pca_plot(all_data, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

# masking based on all data types
sig_list = vector("list", length(all_data))
names(sig_list) = names(all_data)

for(i in 1:length(sig_list)) {
  limma_res = limma_diff(all_data[[i]]$res, groups=list(1:16, 17:50), thresh=5000)
  limma_res$gene = roi$gene[as.numeric(rownames(limma_res))]
  limma_gsea = arrange(limma_res, desc(diff)) %>% filter(adj.P.Val<0.0001, abs(diff)>1) %>% select(gene, diff)
  write_tsv(limma_gsea, paste0("c:/Downloads/tmp/gsea/gene_ranks/otar_", i, ".rnk"), col_names=FALSE)  
  sig_list[[i]] = limma_gsea
}

venn.diagram(lapply(sig_list, function(x) x$gene), filename="out.png", imagetype="png")
venn_res = calculate.overlap(lapply(sig_list, function(x) x$gene))
venn_genes = c(venn_res$a6, venn_res$a7, venn_res$a5, venn_res$a2)
feature_ix = match(venn_genes, roi$gene)

all_data_sliced = all_data
for(i in 1:length(all_data_sliced)) {
  all_data_sliced[[i]]$res = all_data_sliced[[i]]$res[,-feature_ix]
}

pca_data = prep_for_pca_plot(all_data_sliced, annot_1=group_labels, annot_2=single_labels, marks=c("rna","H3K27ac", "H3K4me3", "H3K27me3"))

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
dev.off()

