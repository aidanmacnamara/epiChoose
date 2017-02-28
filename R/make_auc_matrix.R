#' @title TO ADD
#' @description This function takes as input a .csv file with fields: Label, Mark, Bigwig
#' @description This is a new line ...
#' @details What's this?
#' @param input_data TO ADD
#' @param root TO ADD
#' @param roi Region of interest, a Granges object that specifies the genomic regions to calculate AUC
#' @param mark Which mark in \code{input_data} to use

make_auc_matrices <- function(input_data, roi, mark) {
  
  require(tidyverse)
  
  dat = read_csv(input_data)
  
  
}

