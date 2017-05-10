#' @title Make an n*p AUC matrix for a given sample set and mark
#' @description This function takes as input a .csv file with fields: Label, Mark, Bigwig
#' @description This is a new line ...
#' @details What's this?
#' @param input_data The metadata for each sample. A .csv file that contains 'Label', 'Mark', and 'Bigwig' fields.
#' @param root TO ADD
#' @param roi Region of interest, a sorted Granges object that specifies the genomic regions to calculate AUC. It also need 'gene' in the metadata
#' @param mark Which mark in \code{input_data} to use

make_auc_matrix <- function(input_data, roi, mark, tmp_dir, quantile_norm=TRUE) {
  
  require(tidyverse)
  require(stringr)
  require(readr)
  require(limma)
  
  dat = readr::read_csv(input_data)
  dat = dplyr::filter(dat, Mark==mark)
  
  res_matrix = matrix(NA, nrow=dim(dat)[1], ncol=length(roi))
  roi_df = as.data.frame(roi)[,1:3]
  
  colnames(res_matrix) = str_replace_all(apply(roi_df, 1, function(x) paste(as.character(x), collapse="_")), "[[:blank:]]+", "")
  rownames(res_matrix) = dat$Label
  
  tmp_bed = str_replace(tempfile(tmpdir=tmp_dir, pattern="", fileext=".bed"), "\\\\", "")
  roi_df[,2] = as.character(roi_df[,2])
  roi_df[,3] = as.character(roi_df[,3]) # make sure scinot is not used in export
  write_tsv(roi_df, tmp_bed, col_names=FALSE)
  
  for(i in 1:dim(dat)[1]) {
    
    if(is.na(dat$Bigwig[i])) {
      print(paste("Sample", i, "has no data for this mark, skipping ..."))
      next
    } else {
      print(paste("Sample", i))
    }
    
    cmd = paste0("WiggleTools/wiggletools apply AUC ", tmp_bed, " ", dat$Bigwig[i])
    res = system(cmd, intern=TRUE)
    
    if(!length(res)) {
      print(paste("Sample", i, "did not work with WiggleTools - no data returned, skipping ..."))
      next
    }

    res_idx = grep("fixedStep", res)
    
    # need to clean wiggletools output because of, for example:
    # 'fixedStep chrom=chr11 start=867501 step=1' at certain lines
    if(length(res_idx)) res = res[-res_idx]
    # parse res
    res_coords = str_replace(res, "^(.*)\t(.*)\t(.*)\t(.*)$", "\\1_\\2_\\3")
    res_score = as.numeric(str_replace(res, "^(.*)\t(.*)\t(.*)\t(.*)$", "\\4"))
    res_match = match(res_coords, colnames(res_matrix))
    res_score = res_score[res_match]
    res_matrix[i,] = res_score
    
  }
  
  if(quantile_norm) {
    r_ix = which(rowSums(res_matrix, na.rm=TRUE)==0)
    res_trans = t(res_matrix[-r_ix,])
    res_trans = normalizeQuantiles(res_trans)
    res_matrix[-r_ix,] = t(res_trans)
  }
  
  return(list(res=res_matrix, annot=dat))
  
}

