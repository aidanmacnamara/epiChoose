
plot_data_heatmap <- function(mofa_model, plot_view, my_factor, sample_ix, features=50, include_weights=FALSE, transpose=FALSE, ...) {
  
  W <- getWeights(mofa_model)[[plot_view]][,my_factor]
  Z <- getFactors(mofa_model)[,my_factor]
  Z <- Z[sample_ix]
  Z <- Z[!is.na(Z)]
  
  data <- getTrainData(mofa_model, plot_view)[[1]]
  data <- data[,names(Z)]
  data <- data[,apply(data, 2, function(x) !all(is.na(x)))]
  if(is(features, "numeric")) {
    features <- names(tail(sort(abs(W)), n=features))
  }
  data <- data[features,]
  
  order_samples <- names(sort(Z, decreasing=TRUE))
  order_samples <- order_samples[order_samples %in% colnames(data)]
  data <- data[,order_samples]
  
  if(include_weights == TRUE) {
    anno <- data.frame(row.names=names(W[features]), weight=W[features])
    pheatmap(t(data), annotation_row=anno, ...)
  }
  else {
    pheatmap(t(data), ...)
  }
  
}


