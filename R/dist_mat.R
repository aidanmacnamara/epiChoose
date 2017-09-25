#' @title Computes the final distance matrix for comparisons of choice
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param all_data
#' @param comp_ix
#' @return TO ADD

dist_mat <- function(x, comp_ix, labels, plot_labels="", plot_res=TRUE, my_title="", use_corr=FALSE, font_size=15, label_size=4) {
  
  
  get_dists <- function(x) {
    x = x$res
    complete_ix = which(apply(x, 1, function(x) !all(is.na(x)))) # find samples with any data
    
    if(use_corr) {
      y = matrix(NA, nrow=dim(x)[1], ncol=dim(x)[1])
      colnames(y) = rownames(x)
    } else {
      y = matrix(NA, nrow=dim(x)[1], ncol=2)
    }    
    
    rownames(y) = rownames(x)
    x = cor(t(x[complete_ix,]), use="pairwise.complete.obs")
    
    if(use_corr) {
      y[complete_ix,complete_ix] = x
    } else {
      x = 0.5 * (1-x)
      x = cmdscale(x, eig=T)$points
      y[complete_ix,] = x
    }
    
    return(y)
    
  }
  
  all_dists = lapply(x, get_dists)
  
  res = data.frame()
  
  for(k in 1:length(x)) {
    
    out_mat = matrix(NA, nrow=length(comp_ix[[2]]), ncol=length(comp_ix[[1]]))
    for(i in 1:dim(out_mat)[1]) {
      for(j in 1:dim(out_mat)[2]) {
        
        if(use_corr) {
          out_mat[i,j] = 1 - all_dists[[k]][comp_ix[[1]][j],comp_ix[[2]][[i]]]
        } else {
          out_mat[i,j] = sqrt((all_dists[[k]][comp_ix[[1]][j],1] - median(all_dists[[k]][comp_ix[[2]][[i]],1], na.rm=TRUE))^2 + (all_dists[[k]][comp_ix[[1]][j],2] - median(all_dists[[k]][comp_ix[[2]][[i]],2], na.rm=TRUE))^2)
        }
        
      }
    }
    
    res = rbind(res, out_mat)
    
  }
  
  names(res) = labels[comp_ix[[1]]]
  rownames(res) = names(x)
  
  if(plot_res) {
    res_melt = melt(t(res))
    names(res_melt) = c("Cell", "Assay", "Distance")
    label_ix = c()
    for(i in 1:length(plot_labels)) {
      label_ix = c(label_ix, grep(plot_labels[i], res_melt$Cell))      
    }
    res_melt$Cell[-label_ix] = NA
    print(ggplot(res_melt, aes(x=Assay, y=Distance, color=Assay)) + theme_thesis(font_size) + geom_jitter(width=0.1, height=0, shape=17) + geom_text_repel(aes(label=Cell), fontface="bold", size=label_size, force=0.5) + ggtitle(my_title))
  }
  
  return(res)
  
}

