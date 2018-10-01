#' @title Returns the intersect of samples across all ChIP datasets contained in the list
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

prep_across_datatypes <- function(dat_chip) {
 
  # what are the common labels
  all_labels = sort(unique(unlist(lapply(dat_chip, function(x) x$annot$Label))))
  all_ix = lapply(dat_chip, function(x) match(all_labels, x$annot$Label))
  
  dat_chip_filtered = dat_chip
  
  for(i in 1:length(dat_chip)) {
    dat_chip_filtered[[i]]$res = dat_chip[[i]]$res[all_ix[[i]],]
    rownames(dat_chip_filtered[[i]]$res) = all_labels
    dat_chip_filtered[[i]]$annot = dat_chip[[i]]$annot[all_ix[[i]],]
  }
  
  return(dat_chip_filtered)

}


