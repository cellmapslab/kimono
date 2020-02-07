
#' Calculate tau based on frobenius norm
#' 
#' @parameter x matrix/dataframe. must have at least rows & cols > 2
#' @return Frobenius norm of x
#' @examples 
#' x <- matrix(1:20,5,4)
#' calc_tau(x)
calc_tau <- function(x){
  fr_norm <- calc_frobenius_norm(x)
  10^(-1 *fr_norm) 
}

#' Calculate alpha based on frobenius norm and group information
#' 
#' @parameter x matrix/dataframe. must have at least rows & cols > 2
#' @parameter groups, vector of group numbers
#' @return Frobenius norm of x
#' @examples 
#' x <- matrix(1:20,5,4)
#' group <- c(1,1,2,2)
#' calc_alpha(x, group)
calc_alpha <- function(x, group){

  fr_norm <- c()
  for (g in unique(group)) { 
    tmp <-  calc_frobenius_norm(as.matrix(x[, which(group == g)]))
    fr_norm <- c( fr_norm, tmp )
  }
  10^mean( -1 * fr_norm, na.rm = T) 

}

#' Calculate lambda1.se model for crossvalidation
#' 
#' @parameter lambdas vector of lambdas
#' @parameter cv_result, matrix with nrow == folds and ncol ==  length(lamdas)
#' @return lambda1.se
calc_lambda1.se <- function(lambdas,error_cv){
  
  cv <- colMeans(error_cv)
  std         <-  sd(cv) / sqrt(sum(!is.na(cv)))
  lambdas[ min(which(cv <= (min(cv) + std) )) ]
}

#' Extracts groups numbers
#' 
#' @parameter col names vector like 
#' @return vector of numeric groups
#' @examples 
#' names <- c("methylation___cg0123123","rppa___MTOR")
#' infer_prior_groups(names)
parse_prior_groups <- function(names, sep="\\___"){
  as.numeric( as.factor( do.call( rbind, strsplit( names , split = sep ))[,1]))
}

#' detects non informative features 
#' 
#' @parameter matrix x
#' @return vector of bool
#' @export
is_underpowered <- function(x){
  apply( x , 2, var) == 0
}

#' detects non informative features 
#' 
#' @parameter Y_hat predicted matrix with one y_hat per column
#' @return vector of mse
calc_cv_error <- function(Y_hat,y){

  apply( Y_hat , 2, calc_mse, y ) 
}

#' detects non informative features 
#' 
#' @parameter Y_hat predicted matrix with one y_hat per column
#' @return vector of mse
parsing_name <- function(string,sep="___"){
  
  idx <- grepl("___",string) #might be variables like intercept which do not have any prefix
 
  result <- data.frame('prefix'= as.character(string), "id"= as.character(string), stringsAsFactors = FALSE) 
  split_string <- do.call(rbind,strsplit(as.character(string[idx]),'___'))
  result[-idx,1] <- split_string[,1]
  result[-idx,2] <- split_string[,2]
  
  result
}

#' Calculate crossvalidation to estimate lambda1.se & identify potential underpowered features
#' 
#' @parameter y data.table - feature to predict
#' @parameter X data.table - input features with prior names attached to features
#' @parameter model - input features with prior names attached to features
#' @parameter intercept boolean 
#' @parameter seed_cv int - remove randomness for crossvalidation
#' @parameter folds_cv int - defining the amount of folds
#' @return list containing lambda1.se and feature names excluding underpowered ones
calc_cv_sgl <- function(y, x, model = "sparse.grp.lasso", intercept = TRUE, seed_cv = 1234, folds_cv = 5){

  error_cv    <- matrix(NA, ncol=100, nrow = folds_cv) # matrix storing cv errors. oem tests for 100 lambda values
  repeat_cv   <- TRUE
  
  while (repeat_cv) {
    
   
    if(ncol(x) < 2) # in case we exclude all features return empty list
      return(list()) 
    
    set.seed(seed_cv) # set seed to reproduce cv folds
    fold_idx <- sample( rep( 1:folds_cv, length.out = nrow(x) ) )
    
    for (fold in 1:folds_cv) {
      
      #define test and training sets
      test_x  <- x[which(fold_idx == fold), ]    
      test_y  <- y[which(fold_idx == fold)]
      
      train_x <- x[which(fold_idx != fold), ]  
      train_y <- y[which(fold_idx != fold)]
      
      #remove underpowered features by testing the variance in training and test set
      underpowered <- is_underpowered( test_x) |  is_underpowered( train_x)      

      x <- x[ , !underpowered ]
      
      #restart cv if underpowered features were detected
      if( any( underpowered ))
        break   
         
      #transform data table to matrix for calculations
      train_x <- as.matrix(train_x)
      train_y <- as.matrix(train_y)
      test_x  <- as.matrix(test_x)
      test_y  <- as.matrix(test_y)
      
      #prepare oem parameters
      group <- parse_prior_groups(colnames(x))
      tau   <- calc_tau(train_x )
      alpha <- calc_alpha(train_x, group)
      
      # supressing warnings since oem warns for all "p > n -> it might be slow for small sample sizes". 
      # however still the best package for sparse group lasso
      fit_cv <- suppressWarnings( 
        oem(x = train_x ,
            y = train_y , 
            penalty =  model,
            alpha = alpha,
            tau = tau ,
            groups = group, 
            intercept = intercept)
       )

      Y_hat             <- predict( fit_cv, test_x, type = "response" ) # predict test set
      error_cv[fold, ]  <- calc_cv_error(Y_hat,test_y )                 # evaluate test error 
    }
    
    #stop loop if no features got removed without error
    repeat_cv <- ifelse(ncol(x) == ncol(train_x),  FALSE, TRUE) 
  }
  
  lambda1.se  <- calc_lambda1.se(fit_cv$lambda[[1]], error_cv)
  features    <- colnames(x)
  
  list("lambda1.se" = lambda1.se,
       "features" = features )
}

#' Calculate crossvalidation to estimate lambda1.se & identify potential underpowered features
#' 
#' @parameter y data.table - feature to predict
#' @parameter X data.table - input features with prior names attached to features
#' @parameter model string - which model to train. currently only sparse group lasso tested
#' @return edge list for a given input y and x
train_kimono_sgl  <- function(y, x, model = "sparse.grp.lasso", intercept = TRUE, ...){ 
  
  y <- scale(y)
  x <- scale(x)
  
  # estimate best lambda and identify underpowered features 
  cv_result   <- calc_cv_sgl(y, x )
  
  if(length(cv_result)==0) return(c()) #exit function hyere if all features got excluded
  
  # parse cv results
  feature_set <- cv_result$features
  lambda1.se  <- cv_result$lambda
  
  # exclude underpowered features and convert input to matrix
  x <- x[,colnames(x) %in% feature_set, drop = FALSE]
  x <- as.matrix(x)
  y <- as.matrix(y) 
  
  # get sgl input parameters
  group <- parse_prior_groups(colnames(x))
  tau   <- calc_tau(x)
  alpha <- calc_alpha(x, group)
  
  #fit model with supressed warnings on whole dataset. 
  # supressing warnings since oem warns for all "p > n -> it might be slow for small sample sizes". 
  # however still the best package for sparse group lasso
  fit <-suppressWarnings( 
    oem(x = x ,
        y = y , 
        penalty =  model,
        lambda = lambda1.se ,
        alpha = alpha,
        tau = tau ,
        groups = group, 
        intercept = intercept #oem performs better with intercept even though
                              # we scaled the input data. Inferred intercepts
                              # can be ignored since they are almost 0.
        )
  )
  
  y_hat <- predict(fit, x, type = "response")
  
  covariates  <- rownames(fit$beta$sparse.grp.lasso)
  beta        <- as.vector(fit$beta$sparse.grp.lasso)
  performance <- calc_r_square(y, y_hat )
  mse <- calc_mse(y,y_hat)
  
  #return c() if fit is an intercept only models
  if(!any(beta[covariates != "(Intercept)"] != 0)) return(c())

  
  prefix_covariates <- parsing_name(covariates)
  
  data.frame("predictor"=prefix_covariates$id,
                      "value"=beta,
                      "performance"= performance,
                      "mse"=mse,
                      "relation"=prefix_covariates$prefix
                      )
}

#for debugging purpose only 
#DEBUG <- train_kimono_sgl(y,x)
#cat("performance:",DEBUG[1,3],"\n features:",length(which(DEBUG[,2]!= 0)),"of:",length(DEBUG[,2]),"\n","mse:",DEBUG[1,4])
#model = "sparse.grp.lasso"
#intercept = TRUE
#seed_cv = 1234
#folds_cv = 5
