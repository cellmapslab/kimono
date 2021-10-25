#' @rdname kimono estimate R squared based on ESS/TSS
#' @keywords internal
#' @param y vector of double, assumed to be TRUE
#' @param y_hat vector of double, predicted
#' @return R2 value - double
calc_r_square <- function(y, y_hat){
  sum( (y_hat - mean(y))^2 )  / sum( (y - mean(y) )^2 )
}

#' @rdname kimono estimate mse
#' @keywords internal
#' @param y vector of double, assumed to be TRUE
#' @param y_hat vector of double, predicted
#' @return  mse double
calc_mse <- function(y,y_hat){
  mean((y-y_hat)^2)
}

#' @rdname kimono estimate frobenius norm of a matrix
#' @keywords internal
#' @param x vector of double, assumed to be TRUE
#' @return  frobenius norm double
calc_frobenius_norm <- function(x){
  m    <- cor(as.matrix(x))
  sqrt( sum( m[upper.tri(m)]^2) ) / sqrt( (nrow(m)^2 - nrow(m)) /2 )
}

#' @rdname kimono using the prior information to fetsh the right data for X
#' @keywords internal
#' @param node_name fetch this node from data
#' @param prior_network prior network with more node information
#' @param input_data input data
#' @param prior_missing if prior is missing which?
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

  if(length(prior_missing)>0){
    for (layer in prior_missing) {
      features_prior_missing <- c(colnames(input_data[[layer]]),features_prior_missing)
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

#' @rdname kimono remove na, scale
#' @keywords internal
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


#' @rdname kimono remove na, scale
#' @keywords internal
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


#' @rdname kimono check if data is valid
#' @keywords internal
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


#' @rdname kimono Progressbar function
#' @keywords internal
#' @param iterations - node_iterator length
#' @return results as dataframe
combine_results <- function(iterations){
  if(iterations == 1){
    function(...) {
      rbind(...)
    }
  }else{
    pb <- txtProgressBar(min = 1, max = iterations - 1, style = 3)
    count <- 0
    function(...) {
      count <<- count + length(list(...)) - 1
      setTxtProgressBar(pb, count)
      flush.console()
      rbind(...) # this can feed into .combine option of foreach
    }
  }
}


