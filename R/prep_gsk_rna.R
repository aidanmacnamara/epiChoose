#' @title Prepares the GSK RNA data so it is the same shape and format as the ChIP data
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

prep_gsk_rna <- function(r_data_counts, r_metadata, gene_list, chip_labels, rna_labels, quantile_norm=TRUE) {
  
  require(limma)
  require(SummarizedExperiment)
  
  chip_labels = toupper(chip_labels) # make sure everything is in upper cases
  rna_labels = toupper(rna_labels)
  
  # map genes
  mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters='ensembl_gene_id', values=rownames(r_data_counts), mart=mart_1)
  rownames(r_data_counts) = mapping$hgnc_symbol[match(rownames(r_data_counts), mapping$ensembl_gene_id)]
  
  # normalise
  if(quantile_norm) {
    r_ix = which(rowSums(r_data_counts, na.rm=TRUE)==0)
    res_trans = r_data_counts[-r_ix,]
    res_trans = normalizeQuantiles(res_trans)
    r_data_counts[-r_ix,] = res_trans
  }
  
  # transpose
  r_data_counts = t(r_data_counts)
  
  # remove no gene mapping instances
  r_data_counts = r_data_counts[,colnames(r_data_counts)!=""]
  
  # average genes
  r_data_counts = as.matrix( 
    sapply(unique(colnames(r_data_counts)),
           function(col) rowMeans(r_data_counts[,which(colnames(r_data_counts)==col), drop=FALSE]) 
    )
  )
  
  # make sure genes are the same
  r_data_counts = r_data_counts[,match(gene_list, colnames(r_data_counts))]
  
  # TOADD parse samples so epigenetics and rna match
  # create na rows if sample is not present
  match_ix = match(chip_labels, rna_labels)
  ena_names = r_metadata$`Comment[ENA_RUN]`[match_ix]
  ena_ix_data = match(ena_names, rownames(r_data_counts))
  r_data_counts = r_data_counts[ena_ix_data,]
  ena_ix_annot = match(ena_names, r_metadata$`Comment[ENA_RUN]`)
  r_metadata_filt = r_metadata[ena_ix_annot,]
  rownames(r_data_counts) = chip_labels
  
  return(list(res=r_data_counts, annot=r_metadata_filt))
}

