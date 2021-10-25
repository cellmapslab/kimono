#' estimate R squared based on ESS/TSS
#'
#' @parameter y vector of double, assumed to be TRUE
#' @parameter y_hat vector of double, predicted
#' @return R2 value - double
calc_r_square <- function(y, y_hat){
  sum( (y_hat - mean(y))^2 )  / sum( (y - mean(y) )^2 )
}

#' estimate Mean Squared Error
#'
#' @parameter y vector of double, assumed to be TRUE
#' @parameter y_hat vector of double, predicted
#' @return  mse double
calc_mse <- function(y,y_hat){
  mean((y-y_hat)^2)
}

#' estimate frobenius norm of a matrix
#'
#' @parameter x vector of double, assumed to be TRUE
#' @return  frobenius norm double
calc_frobenius_norm <- function(x){
  m    <- cor(as.matrix(x))
  sqrt( sum( m[upper.tri(m)]^2) ) / sqrt( (nrow(m)^2 - nrow(m)) /2 )
}

