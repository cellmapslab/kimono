#' @rdname kimono kNN imputation
#' @keywords internal
#' @param input_data data.table - with data containing NAs to be imputed
#' @param seed integer - random seed
#' @param ... remain arguments for function impute.knn
#' @return data.table with NAs imputed

knn_imputation <- function(input_data, seed, ...){
  imputed_data <- input_data
  message("KNN imputation + KiMONo")
  # Which data is missing?
  # "gene", "methylation", "cnv", "phenotype"
  for (data.table.name in names(input_data)) {
    mis.val.count <- sum(is.na(input_data[[data.table.name]]))
    message(paste0("# of missing values by ", data.table.name, ": ", mis.val.count))
  }
  message("Starting KNN imputation ...")
  # impute missing values on each level (e.g. gene or cnv data.table) separately
  for (data.table.name in names(input_data)) {
    data.table <- input_data[[data.table.name]]
    mis.val.count <- sum(is.na(data.table))
    if (mis.val.count == 0) {
      message(paste0("No missing data. Skip ", data.table.name))
    }
    else{
      message(paste0("Imputing ", data.table.name))
      data <- as.data.frame(data.table)
      # transform
      # originally, knn.impute requires rows columns as samples
      # this is important since default behavior is any row(e.g. gene)
      # with more than rowmax% (default: 50%) missing are imputed using the overall mean per sample (column)
      # additionally, if any column has more than colmax% (default: 80%) missing data, the program halts and reports an error
      data.t <- t(data)
      dataImp <- impute.knn(data.t, rng.seed = seed,...)
      dataImp <- t(dataImp$data)
      message("Finished")
      imputed_data[[data.table.name]] <- as.data.table(dataImp)
    }
  }
  return(imputed_data)
}

#' @rdname kimono Network Guided Multiple Imputation
#' @keywords internal
#' @param input_data data.table - with data containing NAs to be imputed
#' @param prior_network igraph object - network to be used as reference
#' @param corr_threshold numerical - threshold value for correlation
#' @param m numerical - number of multiple imputations to be performed
#' @param ... remain arguments fro function mice
#' @return data.table with m times multiple imputation outputs per matrix in input_data

ngMice <- function(input_data, prior_network, corr_threshold, m, ...){
  nodes <- c()
  for(k in names(input_data)){
    nodes <- c(nodes,paste0(k,'___',colnames(input_data[[k]])))
  }
  new_network <- induced.subgraph(prior_network,V(prior_network)$name[V(prior_network)$name %in% nodes])
  
  isolated = which(degree(new_network)==0)
  new_network = delete.vertices(new_network, isolated)
  
  imputed_data <- list()
  for(layer in names(input_data)){
    
    if(!any(is.na(input_data[[layer]]))){
      pheno <- list()
      for(i in seq(m)){
        pheno[[i]] <- input_data[[layer]]
      }
      
      imputed_data[[layer]] <- pheno
      
      break
    }
    
    node_names <- V(new_network)$name[V(new_network)$layer == layer]
    node_names <- gsub(pattern = paste0(layer,"___"),replacement = '',x = node_names)
    conserve <- which(colnames(input_data[[layer]]) %in% node_names)
    input_data[[layer]] <- input_data[[layer]][,..conserve]
    
    d <- distances(new_network,paste0(layer,"___",colnames(input_data[[layer]])),mode="all")
    d <- d[,rownames(d)]
    to_constrain <- ifelse(d == 0 | d == Inf,0,1)
    colnames(to_constrain) <- rownames(to_constrain) <- gsub(pattern = paste0(layer,"___"),replacement = "",x =  colnames(to_constrain))
    
    diag(to_constrain) <- 0
    
    to_constrain <- to_constrain[which(rowSums(to_constrain) > 0),which(colSums(to_constrain) > 0)]
    data <- as.matrix(input_data[[layer]])
    result <- sapply(colnames(to_constrain),function(i){
      x <- to_constrain[i,]
      to_correlate <- names(x)[x == 1]
      row_name <- i
      if(length(to_correlate) <= 5){return(x)}
      corr_result <- abs(cor(data[,to_correlate],data[,row_name],use="pairwise.complete.obs"))
      names(corr_result) <- to_correlate
      corr_result <- corr_result[order(corr_result,decreasing = TRUE)]
      
      to_flip <- corr_result < corr_threshold
      if(length(to_correlate) - sum(to_flip,na.rm = TRUE) >= 5){
        x[names(corr_result)[6:length(corr_result)]] <- 0   
      }else{
        x[names(corr_result)[to_flip]] <- 0
        
      }
      return(x)
    })
    data <- data[,colnames(result)]
    #mids <- parlmice_wrapper(data = data,m = 5,cluster.seed=1234,n.core = 5,n.imp.core = 1,method="norm", predictorMatrix=result)
    mids <- mice(data = data, predictorMatrix=result, m = m, ...)
    pre_imputed_data <- lapply(seq(m), function(i) complete(mids, action = i))
    imputed_data[[layer]] <- pre_imputed_data
    
  }
  
  return(imputed_data)
}