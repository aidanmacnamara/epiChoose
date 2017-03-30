#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param chip_data A list with a list entry (matrix and annotation) for each mark
#' @param rna_data A list with matrix and annotation
#' @return TO ADD

chip_rna_int <- function(chip_data, rna_data, roi, marks, chip_labels, project_labels, log_data=FALSE) {
  
  dim_1 = dim(rna_data$res)[1]
  dim_2 = dim(rna_data$res)[2]
  if(length(project_labels)==1) project_labels = rep(project_labels, dim_1)
  
  # prepare chip vectors
  chip_vectors = vector("list", length(chip_data))
  for(i in 1:length(chip_data)) {
    if(log_data) {
      chip_vectors[[i]] = log(as.numeric(t(chip_data[[i]]$res)))
    } else {
      chip_vectors[[i]] = as.numeric(t(chip_data[[i]]$res))
    }
  }
  
  if(log_data) {
    rna_vector = log(as.numeric(t(rna_data$res)))
  } else {
    rna_vector = as.numeric(t(rna_data$res))
  }
  
  all_data = data.frame(
    rna_vector,
    chip_vectors,
    factor(rep(roi$gene, dim_1)),
    factor(rep(chip_labels, each=dim_2)),
    factor(rep(project_labels, each=dim_2))
  )
  names(all_data) = c("rna", marks, "gene", "sample", "project")
  
  if(log_data) all_data[all_data==-Inf] = 0 # because of log
  
  return(all_data)
  
}

