#' @title Shrinks a regulatory data matrix to a gene matrix using the gene list provided
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat
#' @return TO ADD
#' @export

convert_reg_matrix <- function(dat, roi, gene_list, reg_window=0, summ_method=c("mean","max","sum")) {
  
  summ_method = match.arg(summ_method)
  
  start(gene_list) = start(gene_list) - reg_window
  end(gene_list) = end(gene_list) + reg_window
  
  my_ol = findOverlaps(gene_list, roi)
  dat_out = matrix(NA, nrow=dim(dat)[1], ncol=length(gene_list))
  if(!is.null(gene_list$hgnc_symbol)) {
    colnames(dat_out) = gene_list$hgnc_symbol
  } else {
    colnames(dat_out) = gene_list$gene
  }
  rownames(dat_out) = rownames(dat)
  
  for(i in 1:dim(dat_out)[2]) {
    if(summ_method=="mean") {
      dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, mean, na.rm=TRUE)
    }
    if(summ_method=="max") {
      dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, max, na.rm=TRUE)
    }
    if(summ_method=="sum") {
      dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, sum, na.rm=TRUE)
    }
  }
  
  dat_out[is.infinite(dat_out)] = NA
  
  return(dat_out)
  
}

