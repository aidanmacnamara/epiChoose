#' @title Takes the set of all experiments and produces a data frame for multi-panel plotting
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @param all_data
#' @return TO ADD
#' @export

prep_for_plot <- function(all_data, annot_1, annot_2, marks, plot_type=c("pca","mds"), which_pcs=c(1,2), roi=NULL) {
  
  plot_data = data.frame()
  plot_type = match.arg(plot_type)
  
  for(i in 1:length(all_data)) {
    
    dat = all_data[[i]]$res
    
    if(all(is.na(dat))) {
      next # there is no data for this data type - skip ...
    }
    
    remove_rows = which(apply(dat, 1, function(x) all(is.na(x)))) # remove samples with no data
    if(length(remove_rows)) {
      dat_na_rm = dat[-remove_rows,]
    } else {
      dat_na_rm = dat
    }
    
    print(paste("Data", i, "has", length(remove_rows), "rows removed."))
    dim(dat_na_rm)
    
    if(toupper(names(all_data)[i])=="RNA") {
      dat_na_rm[is.na(dat_na_rm)] = 0
    }
    
    dat_na_rm = dat_na_rm[,!apply(dat_na_rm, 2, function(x) sd(x)==0)] # remove regions with no variance
    dim(dat_na_rm)
    
    # remove na columns
    dat_na_rm = dat_na_rm[,!apply(is.na(dat_na_rm), 2, all)] # remove regions with no data
    dim(dat_na_rm)
    
    if(is.null(roi)) {
      dat_cut = dat_na_rm
    } else {
      dat_cut = dat_na_rm[,colnames(dat_na_rm) %in% roi]
    }
    dim(dat_cut)
    
    my_coords = matrix(NA, nrow=dim(dat)[1], ncol=2)
    
    if(plot_type=="pca") {
      
      pca_res <- prcomp(dat_cut, scale=TRUE, center=TRUE)
      pca_res_summary <- summary(pca_res)
      
      pc_1 = which_pcs[1]
      pc_2 = which_pcs[2]
      
      if(length(remove_rows)) {
        my_coords[-remove_rows,1] = pca_res$x[,pc_1]
        my_coords[-remove_rows,2] = pca_res$x[,pc_2]
      } else {
        my_coords[,1] = pca_res$x[,pc_1]
        my_coords[,2] = pca_res$x[,pc_2]
      }
      
      to_add = data.frame(my_coords, annot_1, annot_2, pca_res_summary$importance[2,pc_1]*100, pca_res_summary$importance[2,pc_2]*100, marks[i])
      plot_data = rbind(plot_data, to_add)
    }
    
    if(plot_type=="mds") {
      
      d_1 = cor(t(dat_cut))
      d_1 = 0.5 * (1-d_1)
      mds_1 = cmdscale(d_1, eig=T)
      
      if(length(remove_rows)) {
        my_coords[-remove_rows,] = mds_1$points
      } else {
        my_coords = mds_1$points
      }
      
      to_add = data.frame(my_coords, annot_1, annot_2, NA, NA, marks[i])
      plot_data = rbind(plot_data, to_add)
    }
  }
  
  names(plot_data) = c("x", "y", "annot_1", "annot_2", paste0("pc_", which_pcs[1]), paste0("pc_", which_pcs[2]), "mark")
  return(plot_data)
  
}


