#' loading mapping files
#'
#' @param df
#' @param layers - default NA, if no colnames please specify your layer names here
#' @return data frame containing all mappings
load_mapping <- function(df, layers=NA){
  if(length(layers) != 2 ){
    layers <- colnames(df)
  }
  map <- as.data.table(df)
  map <- cbind(map, layers[1], layers[2])
  colnames(map) <- c('A','B','layer_A','layer_B')
  map
}

#' creating a prior network on basis of the prior_map
#'
#' @param prior_map
#' @return igraph network
create_prior_network <- function(prior_map) {

  #cleaning
  idx_rm <- as.character(prior_map$A) %in% c('',NA) |
    as.character(prior_map$B) %in% c('',NA)
  prior_map <- prior_map[!idx_rm,]

  #rm duplicated mapings
  prior_map <- distinct(prior_map)

  #nodes
  A <- distinct(prior_map,A,layer_A)
  colnames(A) <- c('id_name','layer')
  B <- distinct(prior_map,B,layer_B)
  colnames(B) <- c('id_name','layer')
  nodes <- distinct(rbind(A,B))
  nodes <- cbind('id' = paste(nodes$layer,nodes$id_name,sep = '___'),
                 'name' = paste(nodes$layer,nodes$id_name,sep = '___'),
                 nodes)

  #links
  links <- data.table( from = paste(prior_map$layer_A,prior_map$A,sep='___'),
                       to = paste(prior_map$layer_B,prior_map$B,sep='___'),
                       relation = paste(prior_map$layer_A,prior_map$layer_B,sep='___')

  )

  graph_from_data_frame(links, directed = FALSE, vertices = nodes)
}


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
        features_prior_missing <- colnames(input_data[[layer_missing]])
    }
  }

  features <- c(unique(features$name),features_prior_missing)

  for (i in 1:length(input_data)) {
    #extract cross omic relations
    x_idx <- colnames(input_data[[i]]) %in% features
    x <- cbind(x,input_data[[i]][,..x_idx,with = FALSE])
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

  tmp_length <- length(y)

  x <- x[!is.na(y), , drop = FALSE]
  y <- y[!is.na(y), drop = FALSE]

  if(!is.null( dim(x) ) )
    x <- x[ ,!is.na(colSums(x)),drop = FALSE]

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
infer_network <- function(input_data, prior_network,  min_features = 2, sel_iterations = 0, core = 1, specific_layer = NULL, prior_missing = TRUE, DEBUG = FALSE, ...) {


  #preprocess input data
  for (i in names(input_data)) {
    colnames(input_data[[i]]) <- paste(i,colnames(input_data[[i]]),sep='___')
  }


  node_names = V(prior_network)$name
  #check if we only infer a specific layer
  if(!is.null(specific_layer)){
    node_names = V(prior_network)$name[V(prior_network)$layer %in% specific_layer]
  }

  if(DEBUG){
    node_names <- node_names[1:10]
  }

  #TODO: check if number of features are too many for inference
  cl <- parallel::makeCluster(core)
  doParallel::registerDoParallel(cl)
  result <- foreach(node_name = node_names, .combine = 'rbind', .packages = 'kimono')  %dopar% {

    source('~/projects/2020_moni/R_package/git/kimono/R/infer_sgl_model.R')
    source('~/projects/2020_moni/R_package/git/kimono/R/kimono.R')
    source('~/projects/2020_moni/R_package/git/kimono/R/utility_functions.R')

    # can't pass on a node in foreach therefore we have to reselect it here
    #get y and x for a given node
    var_list <- fetch_var(node_name , prior_network, input_data, prior_missing)
    if(is.null(var_list))
      return()

    #remove na and scale data
    var_list <- preprocess_data(var_list$y,var_list$x)

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

    if(possible_error)
      return( )

    if(is.null(subnet))
      return( )

    data.table('target' = parsing_name(node_name)$id , subnet , 'target_layer' = parsing_name(node_name)$prefix)

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
kimono <- function(input_data, prior_network, min_features = 2, sel_iterations = 0 , core = 1, specific_layer = NULL,  ...){

  is_prior_missing <- length(names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]) != 0
  result <- infer_network(input_data, prior_network,  min_features, sel_iterations , core, specific_layer, prior_missing = is_prior_missing, DEBUG =TRUE )

  if( nrow(result) == 0){
    warning('model was not able to infer any associations')
    }else{

      if(is_prior_missing){
        idx_row <- (result$predictor != '(Intercept)' | result[,3] != 0 ) &
          result$predictor_layer %in% names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]
        idx_col <- c('target','predictor','target_layer','predictor_layer')

        tmp <- filter(result,idx_row)[,..idx_col]
        colnames(tmp) <- c('A','B','layer_A','layer_B')

        layer_of_interest <- unique(tmp$layer_B)

        prior_network_new <- create_prior_network(tmp)
        tmp <- infer_network(input_data, prior_network_new,  min_features , sel_iterations , core, specific_layer = layer_of_interest, prior_missing = FALSE )

        result <- rbind(result,tmp)
      }
    }

  result
}


