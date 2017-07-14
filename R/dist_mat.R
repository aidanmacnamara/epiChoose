#' @title Computes the final distance matrix for comparisons of choice
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param all_data
#' @param comp_ix
#' @return TO ADD

dist_mat <- function(x, comp_ix) {
  
  all_dists = lapply(x, function(x) {
    x = x$res
    complete_ix = which(apply(x, 1, function(x) !all(is.na(x))))
    y = matrix(NA, nrow=dim(x)[1], ncol=2)
    rownames(y) = rownames(x)
    x = cor(t(x[complete_ix,]), use="complete.obs")
    x = 0.5 * (1-x)
    x = cmdscale(x, eig=T)$points
    y[complete_ix,] = x
    return(y)
  }
  )
  
  out_all = data.frame()
  
  for(k in 1:length(x)) {
    
    out_mat = matrix(NA, nrow=length(comp_ix[[2]]), ncol=length(comp_ix[[1]]))
    for(i in 1:dim(out_mat)[1]) {
      for(j in 1:dim(out_mat)[2]) {
        out_mat[i,j] = sqrt(
          (
            all_dists[[k]][comp_ix[[1]][j],1] - all_dists[[k]][comp_ix[[2]][i],1]
          )^2 +
            (
              all_dists[[k]][comp_ix[[1]][j],2] - all_dists[[k]][comp_ix[[2]][i],2]
            )^2
        )
      }
    }
    
    out_all = rbind(out_all, out_mat)
    
  }
  
  return(out_all)
  
}


