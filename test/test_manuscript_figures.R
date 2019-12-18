

# FIGURE 1 PCA ------------------------------------------------------------

load("tmp/rld_all.RData")

y = t(rld_all)
dim(y)
y = y[,!apply(y, 2, function(x) all(x==0, na.rm=TRUE))] # remove regions with no variance
dim(y)

pca_res <- prcomp(y, scale=TRUE, center=TRUE)
pca_res_summary = summary(pca_res)
figure_1_pca = data.frame(pca_res$x[,1:2])
names(figure_1_pca) = c("x","y")
figure_1_pca$annot_1 = factor(col_data$cell_type)
figure_1_pca$annot_2 = col_data$condition

save(figure_1_pca, file="~/Dropbox/OTAR020/figures_dat/figure_1_pca.RData") # savepoint
# load("~/Dropbox/OTAR020/figures_dat/figure_1_pca.RData")

ggplot(figure_1_pca, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=4, force=0.5) + theme(legend.position="none")


# FIGURE 1 PILEUP ---------------------------------------------------------

load("tmp/total_signal.RData")
load("tmp/total_diff.RData")
load("tmp/total_gene_orders.RData")
load("tmp/coord_labels.RData")
n_tiles = 21

for(i in 1:length(comps)) {
  
  pdf(paste0("~/Downloads/pile_up_",i,".pdf"), height=10, width=20)
  
  p_1 = filter(total_signal, comp==names(comps)[i]) %>% ggplot(aes(coord_ix, factor(gene, levels=rev(unique(total_gene_orders[[i]]))))) + geom_tile(aes(fill=score_mean)) + facet_wrap(~sample_condition, ncol=length(unique(comps[[i]]$sample_dat$sample_condition))) + theme_thesis(15) + coord_cartesian(xlim=c(2,20)) + xlab("") + scale_fill_gradient(low="white",high="red",limits=c(0,50),na.value="red") + scale_x_continuous(breaks=c(2,11,20), labels=coord_labels[c(2,11,20)]) + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())
  
  p_2 = filter(total_diff, comp==names(comps)[i]) %>% ggplot(aes(x=sample, y=factor(gene, levels=rev(unique(total_gene_orders[[i]]))))) + geom_tile(aes(fill=factor(deregulated))) + theme_thesis(15) + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank()) + coord_cartesian(xlim=c(1.1, length(comps[[i]]$sample_comp)-0.1))
  
  print(plot_grid(p_1, p_2, align="h", rel_widths=c(4,1), axis="bt"))
  
  dev.off()
  
}
