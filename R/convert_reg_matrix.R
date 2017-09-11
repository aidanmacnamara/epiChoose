#' @title Shrinks a regulatory data matrix to a gene matrix using the gene list provided
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat
#' @return TO ADD

convert_reg_matrix <- function(dat, gene_list) {

  my_ol = findOverlaps(gene_list, roi)
  dat_out = matrix(NA, nrow=dim(dat)[1], ncol=length(gene_list))
  colnames(dat_out) = gene_list$hgnc_symbol
  rownames(dat_out) = rownames(dat)
  
  for(i in 1:dim(dat_out)[2]) {
    dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, mean, na.rm=TRUE)
  }
  
  return(dat_out)
  
}

