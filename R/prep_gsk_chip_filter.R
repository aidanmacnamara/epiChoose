#' @title Returns the intersect of samples across all ChIP datasets contained in the list
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

prep_gsk_chip_filter <- function(gsk_chip) {
 
  # what are the common labels
  all_labels = sort(unique(unlist(lapply(gsk_chip, function(x) x$annot$Label))))
  all_ix = lapply(gsk_chip, function(x) match(all_labels, x$annot$Label))
  
  gsk_chip_filtered = gsk_chip
  
  for(i in 1:length(gsk_chip)) {
    gsk_chip_filtered[[i]]$res = gsk_chip[[i]]$res[all_ix[[i]],]
    rownames(gsk_chip_filtered[[i]]$res) = all_labels
    gsk_chip_filtered[[i]]$annot = gsk_chip[[i]]$annot[all_ix[[i]],]
  }
  
  return(gsk_chip_filtered)

}


