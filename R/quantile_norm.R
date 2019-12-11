#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param x Data to normalize (n*p)
#' @return TO ADD

quantile_norm <- function(x) {
  
  r_ix = which(rowSums(x, na.rm=TRUE)==0)
  
  if(is_empty(r_ix)) {
    res_trans = t(x)
    res_trans = normalizeQuantiles(res_trans)
    x = t(res_trans)
  } else {
    res_trans = t(x[-r_ix,])
    res_trans = normalizeQuantiles(res_trans)
    x[-r_ix,] = t(res_trans)
  }
  
  return(x)
  
}

