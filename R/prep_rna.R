#' @title Prepares RNA data so it is the same shape and format as the ChIP data
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param rna_labels
#' @param chip_labels
#' @param col_skip
#' @return TO ADD
#' @export

prep_rna <- function(fpkm_table, gene_list, chip_labels, rna_labels, quantile_norm=TRUE, col_skip=2) {
  
  require(limma)
  require(SummarizedExperiment)
  
  # chip_labels = toupper(chip_labels) # make sure everything is in upper cases
  # rna_labels = toupper(rna_labels)
  
  row_ix = match(chip_labels, rna_labels)
  col_ix = match(gene_list, fpkm_table$`Gene Name`)
  
  fpkm_table = as.matrix(fpkm_table[,-c(1:col_skip)])
  
  rna_output = t(fpkm_table[col_ix,row_ix])
  colnames(rna_output) = gene_list
  rownames(rna_output) = chip_labels
  
  # normalise
  if(quantile_norm) {
    r_ix = which(rowSums(rna_output, na.rm=TRUE)==0)
    if(length(r_ix)) {
      res_trans = t(rna_output[-r_ix,])
      res_trans = normalizeQuantiles(res_trans)
      rna_output[-r_ix,] = t(res_trans)
    } else {
      res_trans = t(rna_output)
      res_trans = normalizeQuantiles(res_trans)
      rna_output = t(res_trans)
    }
  }
  
  return(list(res=data.frame(rna_output, check.names=FALSE), annot=data.frame()))
  
}

