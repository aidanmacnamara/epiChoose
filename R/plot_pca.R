#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

plot_pca <- function(dat, annot_1, annot_2, roi=NULL, which_pcs=c(1,2), out_file) {
  
  require(ggplot2)
  require(ggrepel)
  
  if(length(annot_2)==1) annot_2 = rep(annot_2, length(annot_1))
  
  remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data
  if(length(remove_rows)) dat = dat[-remove_rows,]
  dim(dat)
  
  if(length(remove_rows)) {
    annot_1 = annot_1[-remove_rows]
    annot_2 = annot_2[-remove_rows]
  }
  
  dat = dat[,!apply(dat, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(dat)
  
  # remove na columns
  dat = dat[,!apply(is.na(dat), 2, all)] # remove regions with no data
  dim(dat)
  
  if(is.null(roi)) {
    dat_cut = dat
  } else {
    dat_cut = dat[,colnames(dat) %in% roi]
  }
  dim(dat_cut)
  
  pca_res <- prcomp(dat_cut, scale=TRUE, center=TRUE)
  pca_res_summary <- summary(pca_res)
  
  pc_1 = which_pcs[1]
  pc_2 = which_pcs[2]
  pca_plot <- data.frame(x=pca_res$x[,pc_1], y=pca_res$x[,pc_2], annot_1, annot_2)
  
  png(filename=out_file, height=800, width=1200)
  print(ggplot(pca_plot, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + xlab(paste0("PC", pc_1, ": ", pca_res_summary$importance[2,pc_1]*100, "%")) + ylab(paste0("PC", pc_2, ": ", pca_res_summary$importance[2,pc_2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5))
  dev.off()
  
  return(pca_res)
  
}

