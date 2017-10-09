#' @title Prepares data for a scattplot of data types plotted against each other
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export

spotfire_view <- function(dat, x_axis="H3K27ac", y_axis="RNA", comp_ix=list(c(1),c(3),c(15))) {
  
  dat_1 = data.frame(log(dat[[which(names(dat)==x_axis)]]$res+1))
  dat_2 = data.frame(log(dat[[which(names(dat)==y_axis)]]$res+1))
  
  c_ix = comp_ix[[1]]
  a_ix = comp_ix[[2]]
  t_ix = comp_ix[[3]]
  
  if(length(c_ix)>1) {
    c_dat_1 = apply(dat_1[c_ix,],2,mean,na.rm=TRUE)
    c_dat_2 = apply(dat_2[c_ix,],2,mean,na.rm=TRUE)
  } else {
    c_dat_1 = unlist(dat_1[c_ix,])
    c_dat_2 = unlist(dat_2[c_ix,])
  }
  
  if(length(a_ix)>1) {
    a_dat_1 = apply(dat_1[a_ix,],2,mean,na.rm=TRUE)
    a_dat_2 = apply(dat_2[a_ix,],2,mean,na.rm=TRUE)
  } else {
    a_dat_1 = unlist(dat_1[a_ix,])
    a_dat_2 = unlist(dat_2[a_ix,])
  }
  
  if(length(t_ix)>1) {
    t_dat_1 = apply(dat_1[t_ix,],2,mean,na.rm=TRUE)
    t_dat_2 = apply(dat_2[t_ix,],2,mean,na.rm=TRUE)
  } else {
    t_dat_1 = unlist(dat_1[t_ix,])
    t_dat_2 = unlist(dat_2[t_ix,])
  }
  
  dat_1_out = abs(c_dat_1-a_dat_1) - abs(c_dat_1-t_dat_1)
  dat_2_out = abs(c_dat_2-a_dat_2) - abs(c_dat_2-t_dat_2)
  
  dat_out = data.frame(colnames(dat_1), dat_1_out, dat_2_out)
  names(dat_out) = c("Gene", x_axis, y_axis)
  dat_out = dat_out[apply(dat_out, 1, function(x) !any(is.na(x))),] # remove any points with na
  
  return(dat_out)
  
}


