#' @title Glance at data ...
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export

data_glance <- function(dat) {
  lapply(dat, function(x) x$res[1:5,1:5])
}

