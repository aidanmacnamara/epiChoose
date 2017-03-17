#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

prep_blueprint_rna <- function(quantile_norm=TRUE) {
  
  require(tidyverse)
  require(stringr)
  require(SummarizedExperiment)
  require(biomaRt)
  require(limma)
  
  # http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3827/
  
  load("data/rna/E-MTAB-3827-atlasExperimentSummary.Rdata")
  rna_annot = data.frame(read_tsv("data/rna/E-MTAB-3827.sdrf.txt"), check.names=TRUE)
  rna_dat = experimentSummary[[1]]
  rna_dat_mat = assay(rna_dat)
  
  # map genes
  mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters='ensembl_gene_id', values=rownames(rna_dat_mat), mart=mart_1)
  rownames(rna_dat_mat) = mapping$hgnc_symbol[match(rownames(rna_dat_mat), mapping$ensembl_gene_id)]
  
  # normalise
  if(quantile_norm) {
    r_ix = which(rowSums(rna_dat_mat, na.rm=TRUE)==0)
    res_trans = rna_dat_mat[-r_ix,]
    res_trans = normalizeQuantiles(res_trans)
    rna_dat_mat[-r_ix,] = res_trans
  }
  
  # transpose
  rna_dat_mat = t(rna_dat_mat)
  
  # remove no gene mapping instances
  rna_dat_mat = rna_dat_mat[,colnames(rna_dat_mat)!=""]
  
  # average genes
  rna_dat_mat = as.matrix( 
    sapply(unique(colnames(rna_dat_mat)),
           function(col) rowMeans(rna_dat_mat[,which(colnames(rna_dat_mat)==col), drop=FALSE]) 
    )
  )
  
  # map cell types
  rna_annot_filt = rna_annot[match(rownames(rna_dat_mat), rna_annot$Comment.RUN_NAME.),]
  # rownames(rna_dat_mat) = paste(rna_annot_filt$Comment.donor.ID., rna_annot_filt$Characteristics.cell.type., sep="_") # produces duplicate rownames
  
  return(list(res=rna_dat_mat, annot=rna_annot_filt))
  
}

