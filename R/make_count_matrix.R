#' @title TO ADD
#' @description This function takes as input a .csv file with fields: Label, Mark, Bigwig
#' @description This is a new line ...
#' @details What's this?
#' @param input_data TO ADD
#' @param root TO ADD
#' @param roi Region of interest, a sorted Granges object that specifies the genomic regions to calculate AUC.
#' @param mark Which mark in \code{input_data} to use
#' @return TO ADD
#' @export

make_count_matrix <- function(input_data, roi, mark, change_seqnames=TRUE) {
  
  require(tidyverse)
  require(Rsamtools)
  require(BiocParallel)
  require(DESeq2)
  require(GenomicAlignments) 
  
  if(change_seqnames) {
    e_seqlevels = seqlevels(roi)
    new_seqlevels <- mapSeqlevels(e_seqlevels, "NCBI")
    new_seqlevels <- new_seqlevels[complete.cases(new_seqlevels)]
    roi <- renameSeqlevels(roi, new_seqlevels)
  }
  
  dat = read_csv(input_data)
  dat = dplyr::filter(dat, Mark==mark)
  bamfiles <- BamFileList(dat$Bam, yieldSize=2000000)
  # lapply(bamfiles, seqinfo)
  
  register(MulticoreParam())
  se <- summarizeOverlaps(features=roi, reads=bamfiles, mode="Union", ignore.strand=TRUE)
  
  return(list(res_matrix=t(assays(se)$counts), se=se))
  
}

