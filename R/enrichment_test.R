#' @title Performs enrichment test for input genes against gene sets
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param genes
#' @param genes_ref
#' @return To add
#' @export

enrichment_test <- function(genes, genes_ref, gene_sets, adj="fdr", verbose=FALSE) {
  
  tab = lapply(1:length(gene_sets), function(i) {
    if(verbose==TRUE){
      cat("Processing term", i, names(gene_sets)[i], "\n")
    }
    
    genes_ref = genes_ref[!genes_ref %in% genes]
    
    r_inset = sum(genes_ref %in% gene_sets[[i]])
    r_n_inset = length(genes_ref) - r_inset
    
    g_inset = sum(genes %in% gene_sets[[i]])
    g_n_inset = length(genes) - g_inset
    
    f_mat = matrix(c(g_inset, r_inset, g_n_inset, r_n_inset), nrow=2, ncol=2, byrow=F)
    colnames(f_mat) = c("inSet", "ninSet")
    rownames(f_mat) = c("genes", "reference")
    fish = fisher.test(f_mat, alternative="greater")
    pval = fish$p.value
    inset = r_inset + g_inset
    res = c(g_inset, inset, pval)
    res
  })
  
  r_tab = do.call("rbind", tab)
  r_tab = data.frame(as.vector(names(gene_sets)), r_tab)
  r_tab = r_tab[order(r_tab[, 4]), ]
  padj = p.adjust(r_tab[,4], method=adj)
  tab.out = data.frame(r_tab, padj)
  names(tab.out) = c("Term ID", "Genes", "All", "P Value", "P adjusted")
  
  return(tbl_df(tab.out))
  
}

