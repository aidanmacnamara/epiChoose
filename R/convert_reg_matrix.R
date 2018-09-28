#' @title Shrinks a regulatory data matrix to a gene matrix using the gene list provided
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat
#' @return TO ADD
#' @export

convert_reg_matrix <- function(dat, roi, gene_list, summ_method=c("mean","max","sum","tss","closest"), tss_window=2e3, reg_window=2e3) {
  
  which_max <- function(x) {
    y = which.max(x)
    if(is_empty(y)) {
      return(NA)
    } else {
      return(y)
    }
  }
  
  summ_method = match.arg(summ_method)
  
  if(summ_method=="tss") {
    
    start(gene_list) = gene_list$transcription_start_site - tss_window
    end(gene_list) = gene_list$transcription_start_site + tss_window
    
  } else { # summarise across gene body
    
    start(gene_list) = start(gene_list) - reg_window
    end(gene_list) = end(gene_list) + reg_window
    
  }
  
  my_ol = findOverlaps(gene_list, roi)
  
  dat_out = matrix(NA, nrow=dim(dat)[1], ncol=length(gene_list))
  loc_out = matrix(0, nrow=dim(dat)[1], ncol=dim(dat)[2])
  
  if(!is.null(gene_list$hgnc_symbol)) {
    colnames(dat_out) = gene_list$hgnc_symbol
  } else {
    colnames(dat_out) = gene_list$gene
  }
  rownames(dat_out) = rownames(dat)
  
  for(i in 1:dim(dat_out)[2]) {
    
    loc_ix = subjectHits(my_ol)[queryHits(my_ol)==i] # what are the regulatory indices associated with the gene
    
    if(summ_method=="mean") {
      dat_out[,i] = apply(dat[,loc_ix,drop=FALSE], 1, mean, na.rm=TRUE)
      loc_out[,loc_ix] = 1
    }
    if(summ_method=="max"|summ_method=="tss") {
      dat_out[,i] = apply(dat[,loc_ix,drop=FALSE], 1, max, na.rm=TRUE)
      loc_out[cbind(1:dim(loc_out)[1], loc_ix[unlist(apply(dat[,loc_ix,drop=FALSE], 1, which_max))])] = 1
    }
    if(summ_method=="sum") {
      # dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, function(x) {sum(x, na.rm=TRUE) / sum(width(roi[subjectHits(my_ol)[queryHits(my_ol)==i]]))})
      dat_out[,i] = apply(dat[,subjectHits(my_ol)[queryHits(my_ol)==i], drop=FALSE], 1, sum, na.rm=FALSE)
      loc_out[,loc_ix] = 1
    }
    if(summ_method=="closest") {
      closest_ix = head(order(distance(gene_list[i], roi)), 10) # get the 10 closest regulatory regions around the tss
      dat_out[,i] = apply(dat[,closest_ix,drop=FALSE], 1, sum, na.rm=FALSE)
      loc_out[,closest_ix] = 1
    }
  }
  
  dat_out[is.infinite(dat_out)] = NA
  
  return(list(data=dat_out, locations=loc_out))
  
}

