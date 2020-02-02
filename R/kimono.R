#' extracting the relevant mapping information
#' 
#' @param node , String  
#' @param mapping  ,  data.table with ids in thecolumns 1:2
#' @return vector of feature names
fetch_mappings <- function(node, mapping){
    
  colnames(mapping) <- c('V1','V2')
  #extract mapped features
  features <- rbind(  subset( mapping, V1  == node)[, 2],
                      subset( mapping, V2  == node)[, 1, with=FALSE],
                      use.names=FALSE)
  unique( as.vector( as.matrix( features ) ) )
}

#' using the prior information to fetsh the right data for X
#' 
#' @param node , 
#' @param input_list , 
#' @param mapping_list , 
#' @param main_layer - default = 1
#' @return X sample x feature matrix
fetch_var <- function(node ,input_list, mapping_list, metainfo, main_layer = 1, sep = "___"){  
   
  x <- c()
  
  #iterate over all available mappings to fetch features across all levels 
  for (mapping_idx in 1:length(mapping_list) ) {
    
    mapping_to    <- metainfo$main_to[mapping_idx]  # mapping main layer to X
    mapping_name  <- metainfo$ID[mapping_idx]       # mapping name
  
    all_features <- colnames(input_list[[ mapping_to ]]) # all possible input features
  
    #if a mapping_list is NULL we map all features (i.e. for clinical data)
    if( ncol( mapping_list[[mapping_idx]] ) != 0 ){ 
        features <- fetch_mappings(node , mapping_list[[mapping_idx]] )
    }else{
        features <- all_features 
    }
    
    # if main layer is having a prior to itself it might happen that y is in x
    if( mapping_to == main_layer ) 
      features <- features[features !=  node]    
    
    #check if there are actually features left
    if(length(features) == 0) 
      next
    
    #extract relevant data
    data <- input_list[[ mapping_to ]][, features[features %in% all_features] , with = FALSE ]
    
    #if we have data to add
    if(ncol(data) != 0){
      colnames(data) <- paste0(mapping_name,sep,colnames(data) )
      x <- cbind(x, data)
    }
  }
  
  y <- input_list[[main_layer]][ , node, with=FALSE]
  
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
#' 
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param mapping_list  - list of mappings between each data type one
#' @param metainfo  - table of relation between mappings and input list
#' @param main_layer - which input data set represents the main layer (default = 1)
#' @param min_features - autoexclude models with less than 2 features (default = 2) 
#' @return a network in form of an edge table
infer_network <- function(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2) {

  node_list <- colnames(input_list[[main_layer]]) #iterating ofer node_list of main layer
  
  foreach(node = node_list, .combine = 'rbind')  %do% {
    
    #get y and x for a given node
    var_list <- fetch_var(node, 
                          input_list, 
                          mapping_list, 
                          metainfo)
    
    #remove na and scale data
    var_list <- preprocess_data(var_list$y,var_list$x)
    
    #if not enough features stop here
    if(!is_valid(var_list$x,min_features)) 
      return()
    
    #run model in case the model bugs out catch it     
    possible_error <- tryCatch(
      {
        subnet <- train_kimono_sgl(var_list$y,var_list$x )
        FALSE
      },
      error=function(cond) {TRUE},
      warning=function(cond) {TRUE}
    )
    
    if(possible_error) 
      return( )
    
    if(is.null(subnet)) 
      return( )
    
    data.table('target'=node, subnet)
  }
  
}

#' Run Kimono - Knowledge-guIded Multi-Omic Netowrk inference
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @import dplyr
#' @import foreach 
#' @import oem
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param mapping_list  - list of mappings between each data type one
#' @param main_layer - which input data set represents the main layer (default = 1)
#' @param min_features - autoexclude models with less than 2 features (default = 2) 
#' @param core - if core != 1 kimono will perform an parallell computation
#' @return a network in form of an edge table
#' @export
kimono <- function(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2,...){
  
  result <- infer_network(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2)
  
  if(nrow(result) == 0)
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

