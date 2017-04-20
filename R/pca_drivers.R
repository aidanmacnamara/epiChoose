#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat PCA matrix
#' @return TO ADD

pca_drivers <- function(dat, which_pc=1, thresh=5000) {
  
  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data
  if(length(remove_rows)) {
    dat_na_rm = dat[-remove_rows,]
  } else {
    dat_na_rm = dat
  }
  
  print(paste("Data has", length(remove_rows), "rows removed."))
  dim(dat_na_rm)
  
  dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(dat_na_rm)
  
  # remove na columns
  dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data
  dim(dat_na_rm)
  
  pca_res <- prcomp(dat_na_rm, scale=TRUE, center=TRUE)
  pca_res_summary <- summary(pca_res)
  
  # what regions are driving difference?
  res_l = abs(pca_res$rotation[,1])
  mask_ix = head(order(res_l, decreasing=TRUE), thresh)
  return(mask_ix)
  
}


