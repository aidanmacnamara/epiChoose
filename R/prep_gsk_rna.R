#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD

plot_gsk_rna <- function(labels, gene_list, r_data, r_metadata) {
  
  require(limma)
  
  labels = gsk_chip[[1]]$annot$Label
  gene_list = roi$gene
  load("data/rna/E-MTAB-4101-atlasExperimentSummary.Rdata")
  r_data = experimentSummary
  r_metadata = read_tsv("data/rna/E-MTAB-4101.sdrf.txt")
  
  r_data_counts = assays(r_data[[1]])$counts
  
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
  
  return(list(res=r_data_counts, annot=r_metadata_filt))
}

