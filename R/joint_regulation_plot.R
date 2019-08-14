#' @description TO ADD
#' @details What's this?
#' @return TO ADD
#' @export

joint_regulation_plot <- function(cross_clustering, contrasts=list(list(contrast="c.pma.prime", stimulus="PMA"), list(contrast="c.vd3.prime", stimulus="VD3")), n_sets=3, onlyGeneSets=F, bundle=results_bundle, theme_size=20) {
  
  confusion = padArray(cross_clustering$confusion)
  off_diag = confusion - diag(diag(confusion))
  sets = subset(listTopElements(off_diag)[1:n_sets,], value>0)
  
  if (nrow(sets)<n_sets) {
    warning("Fewer than 'n_sets' diagonals found ...")
  }
  
  sets$cluster = with(sets, paste("(",row,",",col,") ",value, sep=""))
  temp = sets %>% rename(c(col="y", row="x")) %>% inner_join(cross_clustering$allocation, by=c("x","y"))
  
  gene_sets = lapply(split(temp, temp$cluster), function(x) (as.numeric(x$gene)))
  n = as.integer(cross_clustering$n)
  
  g_1 = regulationPlot(contrast=contrasts[[1]]$contrast, stimulus=contrasts[[1]]$stimulus, onlyGeneSets=onlyGeneSets, n=-n, geneSets=gene_sets, bundle=bundle) + labs(title=paste(contrasts[[1]]$contrast, contrasts[[1]]$stimulus, sep="; ")) + theme_thesis(theme_size)
  
  g_2 = regulationPlot(contrast=contrasts[[2]]$contrast, stimulus=contrasts[[2]]$stimulus, onlyGeneSets=onlyGeneSets, n=-n, geneSets=gene_sets, bundle=bundle) + labs(title=paste(contrasts[[2]]$contrast, contrasts[[2]]$stimulus, sep="; ")) + theme_thesis(theme_size)
  
  g_12 = cowplot::plot_grid(g_1, g_2, ncol=2, nrow=1)
  list(geneSets=gene_sets, g_1=g_1, g_2=g_2, g_12=g_12)
  
}

