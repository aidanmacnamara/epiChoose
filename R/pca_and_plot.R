#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param rld DESeq2 rld object
#' @return TO ADD

pca_and_plot <- function(rld, annot_label="none", annot_color="none", annot_shape="none") {
  
  if(class(rld)=="matrix") {
    y = rld
  } else {
    y = t(assays(rld)[[1]])
  }
  
  dim(y)
  y = y[,!apply(y, 2, function(x) sd(x)==0)] # remove regions with no variance
  dim(y)
  
  pca_res <- prcomp(y, scale=TRUE, center=TRUE)
  pca_res_summary = summary(pca_res)
  yy = data.frame(pca_res$x[,1:2])
  names(yy) = c("x","y")
  yy$annot_label = factor(annot_label)
  yy$annot_color = factor(annot_color)
  yy$annot_shape = factor(annot_shape)
  
  my_plot = ggplot(yy, aes(x=x, y=y, color=annot_color)) + geom_point(size=5, aes(shape=annot_shape)) + xlab(paste0("PC", 1, ": ", pca_res_summary$importance[2,1]*100, "%")) + ylab(paste0("PC", 2, ": ", pca_res_summary$importance[2,2]*100, "%")) + theme_thesis() + geom_text_repel(aes(label=annot_label), fontface="bold", size=5, force=0.5) # + theme(legend.position="none")
  
  return(my_plot)
  
}

