#' using the prior information to fetsh the right data for X
#'
#' @param node ,
#' @param input_list
#' @param mapping_list ,
#' @param main_layer - default = 1
#' @return X sample x feature matrix
fetch_var <- function(node_name , prior_network, input_data, prior_missing){

  node <- V(prior_network)[node_name]

  y <- data.table()
  y_idx <- colnames(input_data[[node$layer]]) %in% node_name

  if (any(y_idx)) {
    y <- input_data[[node$layer]][,..y_idx,with = FALSE]
  }else{
    return()
  }

  x <- data.table()
  #identify ID of needed features
  features <- neighbors(prior_network, node$name, mode = c("all"))


  ## get omic layer neighbours
  neighbours_within_layer <- features[features$layer %in% node$layer]
  if( length(neighbours_within_layer) > 0 ){
    for (neighbour in 1:length(neighbours_within_layer)) {
      tmp_features <- neighbors(prior_network, neighbours_within_layer[1], mode = c("all"))
      features <- c(features, tmp_features[!(tmp_features$layer %in% node$layer) ])
    }
  }

  ##check if prior is missing for whole layer
  features_prior_missing <- c()
  if(prior_missing){
    layer_missing <- names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]
    if(length(layer_missing)>0){
      for (layer in layer_missing) {
        features_prior_missing <- c(colnames(input_data[[layer]]),features_prior_missing)
      }
    }
  }

  features <- c(unique(features$name),features_prior_missing)

  for (i in 1:length(input_data)) {
    #extract cross omic relations
    x_idx <- colnames(input_data[[i]]) %in% features
    x <- cbind(x,input_data[[i]][,..x_idx,with = FALSE])
  }

  #remove self loops
  if(any(colnames(y) %in% colnames(x))){
    idx <- which(!colnames(x) %in% colnames(y))
    x <- x[,..idx]
  }


  list("y" = y,
       "x" = x)
}

#' remove na, scale
#'
#' @param y , vector of doubles
#' @param x , matrix features in columns and samples in rows
#' @return x, y without na's
preprocess_data <- function(y, x){

  y <- scale(y)
  x <- scale(x)

  x <- x[which(!is.na(y)), , drop = FALSE]
  y <- y[which(!is.na(y)), drop = FALSE]

  if(!is.null( dim(x) ) )
    x <- x[ ,!is.na(colSums(x)),drop = FALSE]

  list("y"=as.data.table(y),
       "x"=as.data.table(x))
}


#' remove na, scale
#'
#' @param y , vector of doubles
#' @param x , matrix features in columns and samples in rows
#' @return x, y without na's
preprocess_scdata <- function(y, x){

  y <- scale(y)
  x <- scale(x)

  x <- x[which(!is.na(y)), , drop = FALSE]
  y <- y[which(!is.na(y)), drop = FALSE]

  tmp_length <- length(y)


  if(!is.null( dim(x) ) ){
    test <- apply(x, 2,function(y) sum(length(which(is.na(y)))))
    x <- x[,which(colnames(x) %in% names(test[test/tmp_length < 0.3])), drop = FALSE]
    y <- y[which(!is.na(rowSums(x))), drop = FALSE]
    x <- x[which(!is.na(rowSums(x))), , drop = FALSE]
  }


  list("y"=as.data.table(y),
       "x"=as.data.table(x))
}


#' check if data is valid
#'
#' @param min_features , default 5
#' @param x , matrix features in columns and samples in rows
#' @return TRUE / FALSE
is_valid <- function( x, min_features  ){

  if( ncol(x) < min_features )
    return(FALSE)

  if( sum(is.na(colSums( as.matrix(x) ))) > ncol(x)-min_features)
    return(FALSE)

  TRUE
}

#' Infers a model for each node in the main layer
#'
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param prior_network - prior network
#' @param min_features - autoexclude models with less than 2 features (default = 2)
#' @param sel_iterations - run stability selection, defines number of iterations. if !=0 it will perform stability selection
#' @param core - to run network inference in parallel
#' @return a network in form of an edge table
infer_network <- function(input_data, prior_network,  min_features = 2, sel_iterations = 0, core = 1, specific_layer = NULL, prior_missing = TRUE, DEBUG = FALSE, scdata=FALSE, ...) {

  #get all features within the prior network
  node_names <- V(prior_network)$name

  # preprocess data
  # filter out nodes where a prior is missing
  prior_filter <- rep(FALSE,length(node_names))
  for (i in names(input_data)) {
    colnames(input_data[[i]]) <- paste(i,colnames(input_data[[i]]),sep = '___')

    data_filter <- colnames(input_data[[i]]) %in% node_names
    if(any(data_filter))
    {
      input_data[[i]] <- input_data[[i]][,..data_filter]
      #cat(i,'prior coverage ',sum(data_filter)/length(data_filter),'\n')
    }

    prior_filter <- prior_filter | (node_names %in% colnames(input_data[[i]]))
    sum(prior_filter)
  }
  node_names <- node_names[prior_filter]


  #check if we only infer a specific layer
  if(!is.null(specific_layer)){
    node_names <- V(prior_network)$name[V(prior_network)$layer %in% specific_layer]
  }

  if(DEBUG){
    node_names <- node_names[1:100]
  }

  iterations <- length( node_names)

  #TODO: check if number of features are too many for inference
  cl <- parallel::makeCluster(core)
  doParallel::registerDoParallel(cl)
  result <- foreach(node_name = node_names, .combine = 'rbind', .packages = 'kimono')  %dopar% {

    library(igraph)
    library(data.table)
    library(dplyr)
    library(oem)
    library(foreach)
    library(doParallel)
    library(tidyverse)


    # can't pass on an igraph node in foreach therefore we have to reselect it here
    #get y and x for a given node
    var_list <- fetch_var(node_name , prior_network, input_data, prior_missing)
    if(sum(dim(var_list$x))==0)
      return()

    if(sum(dim(var_list$y))==0)
      return()


    #remove na and scale data
    if(scdata){
      var_list <- preprocess_scdata(var_list$y,var_list$x)
    }else{
      var_list <- preprocess_data(var_list$y,var_list$x)
    }

    #if not enough features stop here
    if(!is_valid(var_list$x,min_features))
      return()


    #run model in case the model bugs out catch it
    possible_error <- tryCatch(
      {
        if(sel_iterations != 0){
          subnet <- stability_selection( var_list$y, var_list$x, sel_iterations )
        }
        else{
          subnet <- train_kimono_sgl(var_list$y, var_list$x )
        }
        FALSE
      },
      error = function(cond) {TRUE},
      warning = function(cond) {TRUE}
    )

    if(!exists("subnet"))
      return( )

    if(is.null(subnet))
      return( )

    t <- data.table('target' = parsing_name(node_name)$id , subnet , 'target_layer' = parsing_name(node_name)$prefix)

    return(t)
  }
  parallel::stopCluster(cl)

  return(result)
}


#' Run Kimono - Knowledge-guIded Multi-Omic Netowrk inference
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @import dplyr
#' @import foreach
#' @import oem
#' @import doParallel
#' @import igraph
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param prior_network  - igraph with prior information
#' @param min_features - autoexclude models with less than 2 features (default = 2)
#' @param core - if core != 1 kimono will perform an parallell computation
#' @return a network in form of an edge table
#' @export

kimono <- function(input_data, prior_network, min_features = 2, sel_iterations = 0 , core = 1, specific_layer = NULL, DEBUG = FALSE, scdata=FALSE,  ...){

  time <- Sys.time()
  cat('Run started at : ' , as.character(Sys.time()))
  cat('Input : ')
  for (layers in names(input_data)) {
    cat(layers,': \n')
    cat(' - samples  : ', dim(input_data[[layers]])[1], '\n' )
    cat(' - features : ', dim(input_data[[layers]])[2] ,'\n')
    if( any(layers %in% unique(V(prior_network)$layer)) ){
      cat(' - prior nodes :', sum(V(prior_network)$layer %in% layers))
    }
    cat('\n')
  }

  is_prior_missing <- length(names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]) != 0

  cat('Network inference using prior \n')
  result <- infer_network(input_data, prior_network,  min_features, sel_iterations , core, specific_layer, prior_missing = is_prior_missing, DEBUG, scdata )

  if( nrow(result) == 0){
    warning('KiMONo was not able to infer any associations')
  }else{

    if(is_prior_missing){

      cat('Using inferred effects as priors for ')

      idx_row <- (result$predictor != '(Intercept)' | result[,3] != 0 ) &
                  result$predictor_layer %in% names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]

      idx_col <- c('target','predictor','target_layer','predictor_layer')

      tmp <- filter(result,idx_row)[,..idx_col]
      colnames(tmp) <- c('A','B','layer_A','layer_B')

      layer_of_interest <- unique(tmp$layer_B)

      cat(layer_of_interest, '\n')

      prior_network_new <- create_prior_network(tmp)
      tmp <- infer_network(input_data, prior_network_new,  min_features , sel_iterations , core, specific_layer = layer_of_interest, prior_missing = FALSE )

      result <- rbind(result,tmp)
    }
  }

  cat('Done' , Sys.time() - time)

  result
}


