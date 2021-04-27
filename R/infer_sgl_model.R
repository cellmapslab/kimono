
#' Calculate standard error
#'
#' @parameter x vector
#' @return se of x
calc_se <- function(x, na.rm=FALSE) {
  if (na.rm) 
    x <- na.omit(x)
  sqrt(var(x) / length(x))
} 

#' Calculate tau based on frobenius norm
#'
#' @parameter x matrix/dataframe. must have at least rows & cols > 2
#' @return Frobenius norm of x
calc_tau <- function(x){
  # x <- matrix(1:20,5,4)
  # calc_tau(x)
  fr_norm <- calc_frobenius_norm(x)
  10^(-fr_norm)
}

#' Calculate alpha based on frobenius norm and group information
#'
#' @parameter x matrix/dataframe. must have at least rows & cols > 2
#' @parameter groups, vector of group numbers
#' @return Frobenius norm of x
calc_alpha <- function(x, group){
  # x <- matrix(1:20,5,4)
  # group <- c(1,1,2,2)
  # calc_alpha(x, group)

  fr_norm <- c()
  for (g in unique(group)) {
    tmp <-  calc_frobenius_norm(as.matrix(x[, which(group == g)]))
    fr_norm <- c( fr_norm, tmp )
  }
  10^mean( -fr_norm, na.rm = T)

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
parse_prior_groups <- function(names, sep="\\___"){
  # names <- c("methylation___cg0123123","rppa___MTOR")
  # infer_prior_groups(names)
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

  error_cv    <- matrix(NA, ncol = 100, nrow = folds_cv) # matrix storing cv errors. oem tests for 100 lambda values
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

      x <- x[ , !underpowered, with=FALSE ]

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
#' @parameter intercept - intercept only model
#' @parameter seed - seed of model
#' @return edge list for a given input y and x
train_kimono_sgl  <- function(y, x, model = "sparse.grp.lasso", intercept = TRUE, seed_cv = 1234, ...){

  y <- data.table(scale(y))
  x <- data.table(scale(x))

  # estimate best lambda and identify underpowered features
  cv_result   <- calc_cv_sgl(y, x , seed_cv = seed_cv)

  if(length(cv_result)==0) return(c()) #exit function hyere if all features got excluded

  # parse cv results
  feature_set <- cv_result$features
  lambda1.se  <- cv_result$lambda

  # exclude underpowered features and convert input to matrix
  x <- x[,colnames(x) %in% feature_set, with = FALSE]
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
  r_squared <- calc_r_square(y, y_hat )
  mse <- calc_mse(y,y_hat)

  #return c() if fit is an intercept only models
  if(!any(beta[covariates != "(Intercept)"] != 0)) return(c())


  prefix_covariates <- parsing_name(covariates)

  tibble("predictor"=prefix_covariates$id,
                      "value"=beta,
                      "r_squared"= r_squared,
                      "mse"=mse,
                      "relation"=prefix_covariates$prefix
                      )
}


#' Stability selection and summary of multiple runs
#'
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param nseeds - specifies how many iterations to run
#' @return edge list for a given input y and x, and statistics on multiple runs
stability_selection <- function(y, x, nseeds){

  #Initialize the Seeds and empty Matrix and Vectors for mse and rsquared of length 100
  seeds <- 1:nseeds


  df_values  <- matrix(ncol = nseeds, nrow = ncol(x) + 1 , 0) #+1 accounting for intercept
  df_mse <- rep(NA,nseeds)
  df_r_squared <- rep(NA,nseeds)

  # Groupsparse does not always return all features - therefore we need a vector with all the names
  # and all the relations to calculate the stability selection

  tmp_names <- do.call(rbind , strsplit(colnames(x),'___'))
  idx <- data.frame('id' = 1:(ncol(x) + 1),
                    'merged' = c('(Intercept)___(Intercept)',colnames(x)),
                    'names' = c('(Intercept)',tmp_names[,2]),
                    'relation' = c('(Intercept)',tmp_names[,1])
  )

  # Iterate over each of the seeds and set the seed
  for( i in 1:nseeds){

    # Based on the Seed calculate the Cross Validation to determine the best lambda
    # and calculate the best final fit
    fit <- train_kimono_sgl(y,x, seed_cv = seeds[i])

    # For some seeds groupsparse does not return a model these have to be skipped
    if(is.null(fit)){
      df_mse[,i] <- rep(NA,nseeds)
      df_r_squared[,i] <- rep(NA,nseeds)
      df_values[,i] <- rep(NA,nseeds)

      next
    }

    fit <- cbind(fit, 'merged' = paste(fit$relation,fit$predictor,sep = '___'))
    fit <- merge(idx, fit, by='merged', all.x = TRUE )

    df_mse[i] <-  fit$mse[1]
    df_r_squared[i] <-  fit$r_squared[1]
    df_values[,i] <-  fit$value
  }

  # Dataframe with columns   predictor
  # The value is the frequency of a feature being included in the different seed models
  # Overall R-Squared and MSE are the averages of all R-Squared and MSEs of the different seed models
  # Selected R-Squared and MSE are the averages of the seed models where at least one feature was selected

  tibble(
    'predictor' =  idx$names,
    'relation' =  idx$relation,
    'sel_freq' = (nseeds-apply(df_values,1,function(x){sum(is.na(x))}))/nseeds,
    'mean_value' = rowMeans(df_values, na.rm = T),
    'se_value' = apply(df_values, 1, function(x){calc_se(x, na.rm = T)}),
    'sd_value' = apply(df_values, 1, function(x){sd(x, na.rm = T)}),
    'mean_rsq' = mean(df_r_squared, na.rm = T),
    'sd_rsq' = sd(df_r_squared, na.rm = T),
    'mean_mse' = mean(df_mse, na.rm = T),
    'sd_mse' = sd(df_mse, na.rm = T)
  )
}

#for debugging purpose only
#DEBUG <- train_kimono_sgl(y,x)
#cat("performance:",DEBUG[1,3],"\n features:",length(which(DEBUG[,2]!= 0)),"of:",length(DEBUG[,2]),"\n","mse:",DEBUG[1,4])
#model = "sparse.grp.lasso"
#intercept = TRUE
#seed_cv = 1234
#folds_cv = 5
