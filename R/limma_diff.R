#' @title What are the features that are separating 2 groups?
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat PCA matrix
#' @return TO ADD

limma_diff <- function(dat, groups, thresh=5000) {
  
  require(limma)
  
  dat = t(log(dat))
  dat[is.infinite(dat)] = 0
  
  c_idx = groups[[1]]
  p_idx = groups[[2]]
  colnames(dat) = c(paste("Cell Lines", c_idx, sep="_"), paste("Primary", p_idx, sep="_"))
  samples = data.frame(
    Group=factor(c(rep("Cell Lines", length(c_idx)), rep("Primary", length(p_idx))))
  )
  rownames(samples) = colnames(dat)
  design <- model.matrix(~0 + samples$Group)
  colnames(design) = c("Cell_Lines", "Primary")
  
  fit <- lmFit(dat, design)
  contr <- makeContrasts(diff=Cell_Lines-Primary, levels=design)
  fits <- contrasts.fit(fit, contr)
  ebayes_fits <- eBayes(fits)
  res <- topTableF(ebayes_fits, number=thresh)
  
  # return(as.numeric(rownames(res)))
  return(res)
  
}


