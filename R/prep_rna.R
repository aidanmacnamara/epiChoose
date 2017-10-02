#' @title Prepares RNA data so it is the same shape and format as the ChIP data
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

prep_rna <- function(fpkm_table, gene_list, chip_labels, rna_labels, quantile_norm=TRUE) {
  
  require(limma)
  require(SummarizedExperiment)
  
  # chip_labels = toupper(chip_labels) # make sure everything is in upper cases
  # rna_labels = toupper(rna_labels)
  
  rna_output = matrix(NA, nrow=length(chip_labels), ncol=length(gene_list))
  colnames(rna_output) = gene_list
  rownames(rna_output) = chip_labels
  
  row_ix = match(rna_labels, chip_labels)
  col_ix = match(gene_list, fpkm_table$`Gene Name`)
  
  rna_output[row_ix,] = t(fpkm_table[col_ix,3:dim(fpkm_table)[2]])
  
  # normalise
  if(quantile_norm) {
    r_ix = which(rowSums(rna_output, na.rm=TRUE)==0)
    res_trans = t(rna_output[-r_ix,])
    res_trans = normalizeQuantiles(res_trans)
    rna_output[-r_ix,] = t(res_trans)
  }
  
  return(list(res=data.frame(rna_output, check.names=FALSE), annot=data.frame()))
  
}

