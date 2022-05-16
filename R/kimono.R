#' @rdname kimono Infers a model for each node in the main layer
#' @keywords internal
#' @param input_list - list of omics data. First list element will be used as predictor
#' @param prior_network - prior network
#' @param min_features - autoexclude models with less than 2 features (default = 2)
#' @param sel_iterations - run stability selection, defines number of iterations. if !=0 it will perform stability selection
#' @param core - to run network inference in parallel
#' @param specific_layer - run only on one specific layer
#' @param prior_missing - is prior missing
#' @param scdata - if it is sc data
#' @param saveintermediate - saveintermediate
#' @return a network in form of an edge table
infer_network <- function(input_data, prior_network, imputed_data, method, ADW_calculation, min_features = 2, sel_iterations = 0, core = 1, specific_layer = NULL, prior_missing, scdata=FALSE, saveintermediate = FALSE, ...) {

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

  iterations <- length(node_names)

  cat( as.character(Sys.time()), 'starting inference of ', iterations, ' models\n')
  pb <- txtProgressBar(0, iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  #TODO: check if number of features are too many for inference
  cl <- makeCluster(core)
  registerDoSNOW(cl)
  #result <- foreach(node_name = node_names, .combine = 'rbind', .packages = 'kimono', .options.snow=opts)  %dopar% {
  result <- foreach(node_name = node_names, .combine = 'rbind') %do% {

    subnet <- NULL
    var_list <- NULL

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
      var_list <- preprocess_data(var_list$y,var_list$x, method)
    }

    #if not enough features stop here
    if(!is_valid(var_list$x,min_features,method))
      return()

    #run model in case the model bugs out catch it
    possible_error <- tryCatch(
      {
        if(sel_iterations != 0){
          subnet <- stability_selection( var_list$y, var_list$x, sel_iterations )
        }
        else if(method == 'galasso'){
          
          subnet <- galasso(var_list$y, var_list$x,imputed_data,ADW_calculation)
          
        }else if(method == 'salasso'){
          subnet <- salasso(var_list$y, var_list$x,imputed_data,input_data,ADW_calculation)
        }
        else if (method == "sgl"){
          subnet <- train_kimono_sgl(var_list$y, var_list$x )
        } else {
          #browser()
          subnet <- train_kimono_lasso(x = var_list$x, y = var_list$y,  method = method,...)
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

  if(saveintermediate){
    save_kimono(result)
  }
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
#' @param input_data - list of omics data. First list element will be used as predictor
#' @param prior_network  - igraph with prior information
#' @param min_features - autoexclude models with less than 2 features (defaul
#' @param core - if core != 1 kimono will perform an parallell computation
#' @param specific_layer - specify layer
#' @param DEBUG - only run on subset
#' @param scdata - if single cell data
#' @param infer_missing_prior - if kimono should infer missing prior
#' @param saveintermediate - save each intermediate result
#' @return a network in form of an edge table
#' @export

kimono <- function(input_data, prior_network, imputed_data = NULL, method = 'sgl', ADW_calculation = FALSE, min_features = 2, sel_iterations = 0 , core = 1, specific_layer = NULL, scdata=FALSE, infer_missing_prior = FALSE, saveintermediate = FALSE, ...){
  
  checkmate::assertChoice(method, c("sgl", "lasso_coco", "lasso_hm", "lasso_BDcoco","salasso","galasso"))
  
  time <- Sys.time()
  #cat('run started at : ' , as.character(Sys.time()),'\n')
  cat('1) input data:\nlayer - samples - features - prior features\n')
  for (layers in names(input_data)) {
    cat(layers,' - ',dim(input_data[[layers]])[1],' - ', dim(input_data[[layers]])[2] ,' - ' )
    if( any(layers %in% unique(V(prior_network)$layer)) ){
      cat(sum(V(prior_network)$layer %in% layers), '\n')
    }else{
      cat('0\n')
    }

    if(dim(input_data[[layers]])[2] > 500 & !any(layers %in% unique(V(prior_network)$layer))){
      warning('KiMONo is not recommended to process missing priors for layers with 500+ features')
      return(c())
    }
  }

  layer_prior_missing <- names(input_data)[!(names(input_data) %in% unique(V(prior_network)$layer))]
  layer_prior <- names(input_data)[(names(input_data) %in% unique(V(prior_network)$layer))]

  cat('\n')
  cat('2) inference:\n for layers ',layer_prior,'\n')
  result <- infer_network(input_data, prior_network,  imputed_data, method, ADW_calculation, min_features,sel_iterations , core, specific_layer, prior_missing = layer_prior_missing, scdata, saveintermediate )

  cat('\n')
  if ( nrow(result) == 0) {
    warning('KiMONo was not able to infer any associations')
    return(c())
  }

  if (length(layer_prior_missing) != 0 & infer_missing_prior) {
    cat('\n3) missing prior \n')
    for (layer_of_interest in layer_prior_missing) {
      cat(layer_of_interest,'\n')

      intra_map <- c()
      if(length(names(input_data[[layer_of_interest]])) > 1){
        cat('within layer\n')
        features <- names(input_data[[layer_of_interest]])

        features <- features %>%
                    combn(2)%>%
                    t %>%
                    data.table

        prior_fully_connected <- load_mapping(features,c(layer_of_interest,layer_of_interest)) %>% create_prior_network()
        intra_map <- infer_network(input_data, prior_fully_connected,  imputed_data, method, ADW_calculation, min_features , sel_iterations , core, specific_layer = layer_of_interest, prior_missing = layer_prior_missing, saveintermediate  )
        cat('\n')

        intra_map <- intra_map %>%
                         filter(.[[3]] != 0 ) %>%  #filter edge effect size = 0
                         filter(predictor != '(Intercept)') #filter intercepts and intercept only models


        idx_col <- c('target','predictor','target_layer','predictor_layer')
        intra_map <- intra_map[,..idx_col]
        colnames(intra_map) <- c('A','B','layer_A','layer_B')
      }

      cat('overall layer\n')
      idx_row <- (result$predictor != '(Intercept)' | result[,3] != 0 ) &
                  result$predictor_layer %in% layer_of_interest

      idx_col <- c('target','predictor','target_layer','predictor_layer')

      tmp <- filter(result,idx_row)[,..idx_col]
      colnames(tmp) <- c('A','B','layer_A','layer_B')

      prior_network_new <- create_prior_network(rbind(tmp,intra_map))

      tmp <- infer_network(input_data, prior_network_new,imputed_data, method, ADW_calculation, min_features , sel_iterations , core, specific_layer = layer_of_interest, prior_missing = layer_prior_missing, saveintermediate  )
      cat('\n')
      result <- rbind(result,tmp)
    }
  }


  cat('\n')
  cat('Done' , round((Sys.time() - time)/60,2) , 'min')

  result
}


