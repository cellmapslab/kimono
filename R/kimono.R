


#' using the prior information to fetsh the right data for X
#'
#' @param node ,
#' @param input_list
#' @param mapping_list ,
#' @param main_layer - default = 1
#' @return X sample x feature matrix
fetch_var <- function(node_name , prior_network, input_data){

  node <- V(prior_network)[node_name]

  y <- data.table()
  y_idx <- colnames(input_data[[node$layer]]) %in% node_name

  if(any(y_idx)){
    y <- input_data[[node$layer]][,..y_idx,with=FALSE]
  }else{
    return()
  }

  x <- data.table()
  features <- neighbors(prior_network, node$name, mode = c("all"))
  for (i in 1:length(input_data)) {
    x_idx <- colnames(input_data[[i]]) %in% features$name
    x <- cbind(x,input_data[[i]][,..x_idx,with=FALSE])
  }

  list("y"=y,
       "x"=x)

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
#' @param mapping_list  - list of mappings between each data type one
#' @param metainfo  - table of relation between mappings and input list
#' @param main_layer - which input data set represents the main layer (default = 1)
#' @param min_features - autoexclude models with less than 2 features (default = 2)
#' @param sel_iterations - run stability selection, defines number of iterations. if !=0 it will perform stability selection
#' @param core - to run network inference in parallel
#' @return a network in form of an edge table
infer_network <- function(input_data, prior_network,  min_features = 2, sel_iterations = 0, core = 1, ...  ) {

  #preprocess input data
  for (i in names(input_data)) {
    colnames(input_data[[i]]) <- paste(i,colnames(input_data[[i]]),sep='___')
  }

  #check if mapping exists for all data types
  prior_missing <- !(names(input_data) %in% unique(V(prior_network)$layer))
  if(any(prior_missing)){
    prior_missing <- names(input_data)[prior_missing]
  }

  #TODO: check if number of features are too many for inference

  cl <- parallel::makeCluster(1)
  doParallel::registerDoParallel(cl)
  result <- foreach(node_name = V(prior_network)$name[1:10], .combine = 'rbind', .packages = 'kimono')  %dopar% {


    library(igraph)
    source('~/projects/2020_moni/R_package/git/kimono/R/infer_sgl_model.R')
    source('~/projects/2020_moni/R_package/git/kimono/R/kimono.R')
    source('~/projects/2020_moni/R_package/git/kimono/R/utility_functions.R')


    # can't pass on a node in foreach therefore we have to reselect it here
    #get y and x for a given node
    var_list <- fetch_var(node_name , prior_network, input_data)

    #remove na and scale data
    var_list <- preprocess_data(var_list$y,var_list$x)

    #if not enough features stop here
    if(!is_valid(var_list$x,min_features))
      return()

    #run model in case the model bugs out catch it
    possible_error <- tryCatch(
      {
        if(sel_iterations != 0){
          subnet <- stability_selection( var_list$y, var_list$x,   sel_iterations )
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

    data.table('target' = strsplit(node_name,'___')[[1]][2], subnet)

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
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param mapping_list  - list of mappings between each data type one
#' @param main_layer - which input data set represents the main layer (default = 1)
#' @param min_features - autoexclude models with less than 2 features (default = 2)
#' @param core - if core != 1 kimono will perform an parallell computation
#' @return a network in form of an edge table
#' @export
kimono <- function(input_data, prior_network, min_features = 2, sel_iterations = 0 , core = 1, ...){

  result <- infer_network(input_data, prior_network,  min_features , sel_iterations , core, ...)

  if( nrow(result) == 0)
    warning('model was not able to infer any associations')


  result
}

#source(file = 'R_package/kimono/R/infer_sgl_model.R')
#source(file = 'R_package/kimono/R/utility_functions.R')
#main_layer = 1
#node='A1BG'
#min_features = 2
#cores = 1
#start_time <- Sys.time()
#DEBUG <- c()
#DEBUG <-kimono(input_list, mapping_list,metainfo)
#Sys.time() - start_time
#mean(as.matrix(unique(DEBUG[,4,with=FALSE]))) #0.0266

