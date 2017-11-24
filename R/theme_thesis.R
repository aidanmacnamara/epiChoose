#' @title Plotting parameters
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export

theme_thesis <- function(base_size=30, angle_45=TRUE) {
  
  
  my_theme = theme_bw(base_size=base_size) +
    theme(text=element_text(size=base_size, family = ""), axis.title=element_text(size=rel(1)), axis.text=element_text(size = rel(0.75)), axis.title.y=element_text(vjust = 0.3), axis.title.x=element_text(vjust=0.3), plot.title=element_text(size=rel(1.33), vjust=2), legend.title=element_blank(), legend.key=element_blank())
  
  if(angle_45) {
    my_theme = my_theme + theme(axis.text.x=element_text(angle=45, hjust=1))
  }
  
  return(my_theme)
  
}

