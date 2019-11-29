#' @title workspace size
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export

workspace_size <- function() {
  ws <- sum(sapply(ls(envir=globalenv()), function(x) object.size(get(x))))
  class(ws) <- "object_size"
  paste(ws/1e9, "GB")
}
