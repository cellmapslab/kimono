stability_select <- function(x,y, target, nseeds = 50){

  #Initialize the Seeds and empty Matrix and Vectors for mse and rsquared of length 100
  seeds <- 1:nseeds

  mse_values <- rep(NA, nseeds)
  rsquared_values <- rep(NA, nseeds)

  # Groupsparse does not always return all features - therefore we need a vector with all the names
  # and all the relations to calculate the stability selection

  colnames <- colnames(x)
  names <- sapply(colnames, FUN =function(x) {
    split <- unlist(str_split(string = x, pattern = "___")[[1]][2])
    return(split)
  })
  names <- c("(Intercept)", unname(names))
  art <- sapply(colnames, FUN =function(x) {
    split <- unlist(str_split(string = x, pattern = "___")[[1]][1])
    return(split)
  })
  art <- c("(Intercept)", unname(art))



  coef_matrix <- c()

  # Iterate over each of the seeds and set the seed
  for(i in 1:nseeds){


    # Based on the Seed calculate the Cross Validation to determine the best lambda
    # and calculate the best final fit


    fit <- train_kimono_sgl(y,x, seed_cv = seeds[i])
    #print(fit)

    # For some seeds groupsparse does not return a model these have to be skipped
    if(is.null(fit)){
      mse_values[i] <- NA
      rsquared_values[i]<- NA

      next
    }
    fit <- fit %>% rename(r_squared = performance )

    # generate a boolean vector with true values for the features that have an influence
    imp_features <- fit %>% filter(value != 0) %>% pull(predictor)
    fit_features <- names %in% imp_features


    # Calculate the logical vector which features are not 0 in the model
    # bind this vector to the coefficient matrix

    #coef_matrix <- base::cbind(coef_matrix, unname(as.matrix(fit[,"value"]!=0)) )
    coef_matrix <- base::cbind(coef_matrix, fit_features)

    # Calculate R-Squared and MSE by comparing real values to the values predicted by the model
    #y_hat <- as.vector(predict(fit, s =best_lambda, newx=x))
    mse_values[i] <- fit$mse[1]
    rsquared_values[i]<- fit$r_squared[1]

  }

  rownames(coef_matrix) <- names

  #print(coef_matrix)


  # Dataframe with columns of the target and predictor
  # The value is the frequency of a feature being included in the different seed models
  # Overall R-Squared and MSE are the averages of all R-Squared and MSEs of the different seed models
  # Selected R-Squared and MSE are the averages of the seed models where at least one feature was selected

  fit_df <- tibble(
    target = rep(target, nrow(coef_matrix)),
    predictor = rownames(coef_matrix),
    value = rowSums(coef_matrix)/ncol(coef_matrix),
    nr_of_col = ncol(coef_matrix),
    overall_rsq = mean(rsquared_values, na.rm = T),
    selected_rsq = mean(rsquared_values[colSums(coef_matrix) != 1],  na.rm=T),
    overall_mse = mean(mse_values, na.rm= T),
    selected_mse = mean(mse_values[colSums(coef_matrix) != 1], na.rm =T),
  )

  fit_df$relation = art


  return(fit_df)

}






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
#' @param input_list
#' @param mapping_list ,
#' @param main_layer - default = 1
#' @return X sample x feature matrix
fetch_var <- function(node ,input_list, mapping_list, metainfo, main_layer = 1, sep = "___"){

  x <- c()
  #print(paste("Node", node, sep = " "))
  #print(paste("Main Layer ",main_layer))

  #iterate over all available mappings to fetch features across all levels
  for (mapping_idx in 1:length(mapping_list) ) {
    #print(mapping_idx)

    mapping_to    <- metainfo$main_to[mapping_idx]  # mapping main layer to X
    mapping_name  <- metainfo$ID[mapping_idx]       # mapping name

    #print(mapping_to)
    #print(mapping_name)

    #print(input_list[[mapping_to]])

    all_features <- colnames(input_list[[ mapping_to ]]) # all possible input features
    #print(all_features)

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

    #print(features)
    #extract relevant data
    data <- input_list[[ mapping_to ]][, features[features %in% all_features] , with = FALSE ]



    #if we have data to add
    if(ncol(data) != 0){
      colnames(data) <- paste0(mapping_name,sep,colnames(data) )
      x <- cbind(x, data)
    }
  }

  #print(x)

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
infer_network <- function(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2, stab_sel = FALSE) {

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
        if(stab_sel == FALSE){
          subnet <- train_kimono_sgl(var_list$y,var_list$x )
        }
        else{
        subnet <- stability_select(x = var_list$x, y = var_list$y, target = node )
        }
        FALSE
      },
      error=function(cond) {TRUE},
      warning=function(cond) {TRUE}
    )

    if(possible_error)
      return( )

    if(is.null(subnet))
       return( )


    if(stab_sel==F){
    data.table('target'=node, subnet)
    }
    else{
    subnet
    }
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
kimono <- function(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2,stab_sel = FALSE ,...){

  result <- infer_network(input_list, mapping_list, metainfo,  main_layer = 1, min_features = 2, stab_sel = stab_sel)

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

