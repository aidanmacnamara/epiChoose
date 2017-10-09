#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param dat
#' @return TO ADD
#' @export

data_shape <- function(x) {
  
  return(lapply(x, function(y) dim(y[[1]])))
  
}

