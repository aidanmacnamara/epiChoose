#' @title Computes the final distance matrix for comparisons of choice
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param x The original data object, sliced and filtered as necessary
#' @param comp_ix A list of length 2 giving the 1. target(s) and 2. the candidate(s)
#' @return A list containing the ggplot objects and statistics
#' @export

# tmp = dat
# genes = msig_go_bp[[3101]]
# col_ix = which(colnames(tmp[[1]]$res) %in% genes)
# c_cells = rownames(tmp[[1]]$res)[c(205,199,217,211)]
# t_cells = c("S001S7_inflammatory macrophage","S001MJ_inflammatory macrophage","S0022I_inflammatory macrophage")
# c_ix = match(c_cells, rownames(dat[[1]]$res))
# t_ix = match(t_cells, rownames(dat[[1]]$res))
# for(j in 1:length(tmp)) { # each data type
#   tmp[[j]]$res = tmp[[j]]$res[c(c_ix,t_ix),col_ix,drop=FALSE]
# }
# 
# c_ix = match(c_cells, rownames(tmp[[1]]$res))
# t_ix = match(t_cells, rownames(tmp[[1]]$res))
# comp_ix = list(c_ix, t_ix)

dist_mat <- function(tmp, comp_ix, labels, my_title="", font_size=15, label_size=4, label_points=TRUE, get_p_values=FALSE) {
  
  tmp_copy = tmp
  
  if(length(comp_ix[[2]])>1) { # get the mean of the target cells
    for(i in 1:length(tmp_copy)) {
      tmp_copy[[i]]$res = rbind(tmp_copy[[i]]$res, apply(tmp_copy[[i]]$res[comp_ix[[2]],], 2, mean, na.rm=TRUE))
    }
  } else {
    for(i in 1:length(tmp_copy)) {
      tmp_copy[[i]]$res = rbind(tmp_copy[[i]]$res, tmp_copy[[i]]$res[comp_ix[[2]],])
    }
  }
  
  comp_target = dim(tmp_copy[[1]]$res)[1]
  
  get_dists <- function(x) {
    
    x = x$res
    complete_ix = which(apply(x, 1, function(x) !all(is.na(x)))) # find samples with any data
    
    y = matrix(NA, nrow=dim(x)[1], ncol=dim(x)[1])
    colnames(y) = rownames(x)
    rownames(y) = rownames(x)
    
    if(!is_empty(complete_ix)) {
      x = cor(t(x[complete_ix,]), use="pairwise.complete.obs")
      y[complete_ix,complete_ix] = x
    }
    
    return(y)
    
  }
  
  all_dists = lapply(tmp_copy, get_dists)
  res = do.call("rbind", lapply(all_dists, function(x) 1-x[comp_ix[[1]], comp_target]))
  res_melt = melt(t(res))
  names(res_melt) = c("Cell", "Assay", "Distance")
  
  
  get_corr_p <- function(x) {
    x = x$res
    y = rep(NA, length(comp_ix[[1]]))
    if(all(is.na(x[comp_target,]))) {
      return(y)
    } else {
      for(i in 1:length(y)) {
        y[i] = cor.test(x[comp_ix[[1]][i],], x[comp_target,])$p.value
      }
      return(y)
    }
  }
  
  if(get_p_values) {
    corr_p = data.frame(do.call("rbind", lapply(tmp_copy, get_corr_p)))
    names(corr_p) = rownames(tmp[[1]]$res)[comp_ix[[1]]]
  }
  
  
  if(label_points) {
    p_1 = ggplot(res_melt, aes(x=Assay, y=Distance, color=Assay)) + theme_thesis(font_size) + geom_jitter(width=0.1, height=0, shape=17, size=5) + geom_text_repel(aes(label=Cell), fontface="bold", size=label_size, force=0.5) + ggtitle(my_title)
  } else {
    p_1 = ggplot(res_melt, aes(x=Assay, y=Distance, color=Assay)) + theme_thesis(font_size) + geom_jitter(width=0.1, height=0, shape=17, size=5) + ggtitle(my_title)
  }
  
  y_all = data.frame()
  for(j in 1:length(tmp)) {
    x = tmp[[j]]$res
    if(all(is.na(x))) next
    y = melt(as.matrix(x))
    y = cbind(y, names(tmp[j]))
    names(y) = c("Cell Line", "Gene", "Score", "Assay")
    # if(length(which((filter(y, !is.na(Score)) %>% dplyr::select(Gene) %>% table()) == length(unique(y$`Cell Line`)))) < 2) next # if there are < 2 genes with observations from all cell lines, skip ...
    y_all = rbind(y_all, y)
  }
  y_all = tbl_df(y_all)
  
  p_2 = ggplot(y_all, aes(x=Gene, y=Score)) + geom_bar(aes(fill=`Cell Line`), position="dodge", stat="identity") + facet_wrap(~Assay, nrow=2, scales="free")
  
  p_3 = ggplot(y_all, aes(x=`Cell Line`, y=log(Score+1))) + geom_boxplot() + facet_wrap(~Assay, nrow=2, scales="free")
  
  p_4 = ggplot(y_all, aes(x=`Cell Line`, y=Score)) + geom_point() + geom_line(aes(group=Gene)) + facet_wrap(~Assay, nrow=2, scales="free")
  
  if(get_p_values) {
    return(list(plots=list(p_1,p_2,p_3,p_4), stats=list(corr=res, mean=y_all, corr_p=corr_p)))
  } else {
    return(list(plots=list(p_1,p_2,p_3,p_4), stats=list(corr=res, mean=y_all, corr_p=NA)))
  }
  
}

