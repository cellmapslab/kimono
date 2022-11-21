#' @rdname kimono name parsing
#'
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param model string - which model to train. currently only sparse group lasso tested
#' @param folds_cv nr of folds for cv
#' @param nlambdas nr of lambda paramters to be tested
#' @param selection either "lambda.min.index" or "lambda.1se.index"
#'
#' @return edge list for a given input y and x

train_kimono_lasso <- function(x, y, method, folds_cv = 5, seed_cv=1234, nlambdas= 50, selection= "lambda.min.index", rm_underpowered = FALSE){
  set.seed(seed_cv)
  
  x <- x[which(!is.na(y)), , drop = FALSE]
  y <- y[which(!is.na(y)), drop = FALSE]
  
  y <- scale(y)
  x <- scale(as.matrix(x))
  
  if(ncol(x) < 3) # in case we exclude all features return empty list
    return(c())
  
  if(rm_underpowered){
    x <- rm_underpowered_feat(x, folds_cv=folds_cv)
  }
  
  fold_idx <- sample(rep( 1:folds_cv, length.out = nrow(x)))
  #cv_fit <- run_hm(x,y, nlambdas=nlambdas, fold_idx=fold_idx)
  
  if(method == "lasso_coco") {
    cv_fit <- run_coco(x,y, nlambdas=nlambdas, fold_idx=fold_idx)
  } 
  if (method == "lasso_hm") {
    cv_fit <- run_hm(x,y, nlambdas=nlambdas, fold_idx=fold_idx)
  }
  if (method == "lasso_BDcoco"){
    #browser()
    res_coco <- run_BDcoco(x,y,nlambdas=nlambdas)
    cv_fit <- res_coco$cv_fit
    x <-   res_coco$xnew
  }
  
  if(method == "lasso_BDcoco"){
    
    beta <- cv_fit$beta.opt
    
    if(!any(beta != 0)){ return(c())}
    
    y_hat <- x%*% beta
    
    mse <- calc_mse(y,y_hat)
    
    r_squared <- calc_r_square(y,y_hat)
    
    covariates <- cv_fit$vnames
    
    value <- beta
    
  } else {
    
    beta <- cv_fit$fit$beta[,cv_fit[selection][[1]]  ]
    if(!any(beta != 0)){ return(c())}
    
    y_hat <- predict(cv_fit$fit, x)[, cv_fit[selection][[1]] ]
    
    mse <- calc_mse(y,y_hat)
    
    intercept <- cv_fit$fit$a0[, cv_fit[selection][[1]] ]
    
    r_squared <- calc_r_square(y, y_hat )
    
    covariates  <- c("(Intercept)", rownames(cv_fit$fit$beta))
    value <- c(intercept, beta)
  }
  
  prefix_covariates <- parsing_name(covariates)
  
  tibble("predictor" = prefix_covariates$id,
         "value" = value,
         "r_squared" = r_squared,
         "mse" = mse,
         "predictor_layer" = prefix_covariates$prefix,
  )
}

#' @rdname kimono name parsing
#'
#' @param x
#' @param y
#' @param nlambdas
#' @param fold_idx
#'
#' @return

run_BDcoco <- function(x,y, nlambdas){
  
  phenotype_cols <- grep("phenotype___", colnames(x))
  
  if(length(phenotype_cols) < 1) return(list())
  
  x <- cbind(x[,phenotype_cols, drop=FALSE], x[,-phenotype_cols, drop= FALSE])
  
  nr_uncorrupted <- length(phenotype_cols) + 1 # adding one for intercept
  nr_corrupted <- dim(x)[2]-length(phenotype_cols)
  
  x <- cbind(rep(1,nrow(x)),x)
  colnames(x)[1]<- "(Intercept)"
  
  k <- NULL
  
  for(i in c(10:3)){
    if(dim(x)[1]%%i == 0 ) {
      k <- i
      break()
    }
  }
  
  if(is.null(k)) stop("Bug of BDCoCo: sample number needs to have a devivder between 10 and 3")
  
  cv_fit <- BDcocolasso::coco(Z = x, y = y, n=dim(x)[1], p=dim(x)[2], p1=nr_uncorrupted, p2=nr_corrupted,
                              step = nlambdas, K=k,tau=NULL, etol = 1e-4, mu = 10, center.y = FALSE,
                              noise="missing", block= TRUE, penalty= "lasso", mode = "ADMM")
  
  return(list(cv_fit=cv_fit, xnew= x))
}

#' @rdname kimono cocolasso runner
#'
#' @param x
#' @param y
#' @param nlambdas
#' @param fold_idx
#'
#' @return
#' @export
#'
#' @examples

run_coco <-  function(x,y, nlambdas, fold_idx){
  cv_fit <- hmlasso::cv.hmlasso(x, y,nlambda=50, seed = 1234, lambda.min.ratio=1e-1,
                                foldid=fold_idx, direct_prediction=TRUE,
                                positify="admm_max", weight_power = 0)
  return(cv_fit)
}

#' @rdname kimono hmlasso runner
#'
#' @param x
#' @param y
#' @param nlambdas
#' @param fold_idx
#'
#' @return

run_hm <- function(x,y, nlambdas, fold_idx){
  
  weight_powers <- c(0.5, 1, 1.5, 2)
  fits <- vector("list", length= length(weight_powers))
  MSEs <- vector("numeric", length= length(weight_powers))
  
  for(i in seq_along(weight_powers)){
    fits[[i]] <- hmlasso::cv.hmlasso(x, y, nlambda=nlambdas, seed = 1234, lambda.min.ratio=1e-1,
                                     foldid=fold_idx, direct_prediction=TRUE,
                                     positify="admm_frob", weight_power = weight_powers[i])
  }
  
  y_hat <- predict(fits[[i]]$fit, x)[, fits[[i]]$lambda.min.index]
  MSEs[i]<- calc_mse(y,y_hat)
  
  return(fits[[which.min(MSEs)]])
}

#' @rdname kimono Stack adaptive lasso
#' @keywords internal
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param imputed_data - imputed data
#' @param input_data - original input data with NAs
#' @param ADW_calculation - if adaptive weight must be calculated
#' @param nseeds - specifies how many iterations to run
#' @return edge list for a given input y and x, and statistics on multiple runs

salasso <- function(var_y,var_x,imputed_data,input_data,ADW_calculation,set_seed=1234){
  #save(var_y,var_x,imputed_data,input_data,ADW_calculation,file='saenet_chechek.RData')
  set.seed(set_seed)
  
  temp <- matrix()
  for(layer in names(input_data)){
    dat <- input_data[[layer]]
    temp <- cbind(temp,dat)
    temp <- temp[,-1]
  }
  
  all_imputed <- list()
  
  for(mat in seq(length(imputed_data[[1]]))){
    tmp <- matrix()
    for(layer in names(imputed_data)){
      dat <- imputed_data[[layer]][[mat]]
      colnames(dat) <- paste0(layer,'___',colnames(dat))
      tmp <- cbind(tmp,dat)
    }
    all_imputed[[mat]] <- tmp[,-1]
  }
  
  comm <- colnames(temp) %in% colnames(all_imputed[[1]])
  temp <- temp[,comm,with=FALSE]
  
  x <- list()
  y <- list()
  
  for(dat in seq(length(all_imputed))){
    pos_x <- which(colnames(all_imputed[[dat]]) %in% colnames(var_x))
    x[[dat]] <- as.matrix(as.data.table(scale(all_imputed[[dat]][,pos_x])))
    pos_y <- which(colnames(all_imputed[[dat]]) %in% colnames(var_y))
    y[[dat]] <- as.vector(scale(all_imputed[[dat]][,pos_y]))
  }
  
  pf <- rep(1, ncol(var_x))
  adWeight <- rep(1, ncol(var_x))
  weights <- 1 - rowMeans(is.na(temp[,colnames(temp) %in% c(colnames(x[[1]]),colnames(y[[1]])), with=FALSE]))
  fit <- cv.saenet(x, y, pf, adWeight,weights = weights,nlambda=50,alpha=1)
  
  if(ADW_calculation){
    message("calculation adaptive weights")
    cv_coef <- coef(fit)
    
    js <- rep(0,ncol(var_x))
    js <- cv_coef[-1]**2 + js
    
    v = log(ncol(var_x))/log(nrow(var_x) * length(all_imputed))
    
    gamma = ceiling(2 * v/1 - v) + 1
    
    ai = (abs(cv_coef[-1]) + 1 / nrow(var_x) * length(all_imputed))**(-gamma)
    
    fit <- cv.saenet(x, y, pf, adWeight=ai,weights = weights,nlambda=50,alpha=1)
    
  }
  
  cv_coef <- coef(fit)
  
  rsquares <- c()
  mses <- c()
  
  for(k in seq(length(all_imputed))){
    y_hat <- apply(x[[k]],1,function(int){ sum(int * cv_coef[-1])})
    y_hat <- y_hat + cv_coef[1]
    
    rsquares <- c(rsquares, calc_r_square(y[[k]],y_hat))
    mses <- c(mses, calc_mse(y[[k]],y_hat))
  }
  
  r2 <- mean(rsquares)
  mse <- mean(mses)
  
  tibble("predictor" = c('(Intercept)',unlist(strsplit(x = colnames(var_x),split = "___"))[seq(2,ncol(var_x)*2,by=2)]),
         "value" = cv_coef,
         "r_squared" = rep(r2,length(cv_coef)),
         "mse" = rep(mse,length(cv_coef)),
         "predictor_layer" = c('(Intercept)',unlist(strsplit(x = colnames(var_x),split = "___"))[seq(1,ncol(var_x)*2,by=2)])
  )
  
}


#' @rdname kimono Stack adaptive lasso
#' @keywords internal
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param imputed_data - imputed data
#' @param input_data - original input data with NAs
#' @param ADW_calculation - if adaptive weight must be calculated
#' @param nseeds - specifies how many iterations to run
#' @return edge list for a given input y and x, and statistics on multiple runs

galasso <- function(var_y,var_x,imputed_data,ADW_calculation,set_seed=1234){
  #save(var_y,var_x,imputed_data,ADW_calculation,file='galasso_chechek.RData')
  set.seed(set_seed)
  
  all_imputed <- list()
  for(mat in seq(length(imputed_data[[1]]))){
    tmp <- matrix()
    for(layer in names(imputed_data)){
      dat <- imputed_data[[layer]][[mat]]
      colnames(dat) <- paste0(layer,'___',colnames(dat))
      tmp <- cbind(tmp,dat)
    }
    all_imputed[[mat]] <- tmp[,-1]
  }
  
  x <- list()
  y <- list()
  
  for(dat in seq(length(all_imputed))){
    x[[dat]] <- as.matrix(as.data.table(scale(all_imputed[[dat]][,colnames(var_x)])))
    y[[dat]] <- as.vector(scale(all_imputed[[dat]][,colnames(var_y)]))
  }
  
  pf <- rep(1, ncol(var_x))
  adWeight <- rep(1, ncol(var_x))
  fit <- cv.galasso(x, y, pf, adWeight,nlambda=50)
  
  if(ADW_calculation){
    message("calculation adaptive weights")
    cv_coef <- coef(fit)
    
    js <- rep(0,ncol(var_x))
    js <- cv_coef[-1]**2 + js
    
    v = log(ncol(var_x) * length(all_imputed))/log(nrow(var_x) * length(all_imputed))
    
    gamma = ceiling(2 * v/1 - v) + 1
    
    ai = (sqrt(js) + 1 / nrow(var_x) * length(all_imputed))**(-gamma)
    
    fit <- cv.galasso(x, y, pf, adWeight=ai,nlambda=50)
    
  }
  
  cv_coef <- coef(fit)
  
  beta_mean <- rep(0,ncol(var_x)+1)
  beta_mean <- cv_coef + beta_mean
  
  
  rsquares <- c()
  mses <- c()
  
  for(k in seq(length(all_imputed))){
    y_hat <- apply(x[[k]],1,function(int){ sum(int * cv_coef[-1])})
    y_hat <- y_hat + cv_coef[1]
    
    rsquares <- c(rsquares, calc_r_square(y[[k]],y_hat))
    mses <- c(mses, calc_mse(y[[k]],y_hat))
  }
  
  r2 <- mean(rsquares)
  mse <- mean(mses)
  
  tibble("predictor" = c('(Intercept)',unlist(strsplit(x = colnames(var_x),split = "___"))[seq(2,ncol(var_x)*2,by=2)]),
         "value" = beta_mean,
         "r_squared" = rep(r2,length(beta_mean)),
         "mse" = rep(mse,length(beta_mean)),
         "predictor_layer" = c('(Intercept)',unlist(strsplit(x = colnames(var_x),split = "___"))[seq(1,ncol(var_x)*2,by=2)])
  )
  
}

#' @rdname kimono  Calculate tau based on frobenius norm
#' @keywords internal
#' @param x matrix/dataframe. must have at least rows & cols > 2
#' @return Frobenius norm of x
calc_tau <- function(x){
  # x <- matrix(1:20,5,4)
  # calc_tau(x)
  fr_norm <- calc_frobenius_norm(x)
  10^(-fr_norm)
}

#' @rdname kimono  Calculate alpha based on frobenius norm and group information
#' @keywords internal
#' @param x matrix/dataframe. must have at least rows & cols > 2
#' @param groups, vector of group numbers
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

#' @rdname kimono  Calculate lambda1.se model for crossvalidation
#' @keywords internal
#' @param lambdas vector of lambdas
#' @param cv_result, matrix with nrow == folds and ncol ==  length(lamdas)
#' @return lambda1.se
calc_lambda1.se <- function(lambdas,error_cv){

  cv <- colMeans(error_cv)
  std         <-  sd(cv) / sqrt(sum(!is.na(cv)))
  lambdas[ min(which(cv <= (min(cv) + std) )) ]
}

#' @rdname kimono  Extracts groups numbers
#' @keywords internal
#' @param col names vector like
#' @return vector of numeric groups
parse_prior_groups <- function(names, sep="\\___"){
  # names <- c("methylation___cg0123123","rppa___MTOR")
  # parse_prior_groups(names)
  as.numeric( as.factor( do.call( rbind, strsplit( names , split = sep ))[,1]))
}

#' @rdname kimono  detects non informative features
#' @keywords internal
#' @param matrix x
#' @return vector of bool
#' @export
is_underpowered <- function(x){
  apply( x , 2, var) == 0
}

#' @rdname kimono  detects non informative features
#' @keywords internal
#' @param Y_hat predicted matrix with one y_hat per column
#' @return vector of mse
calc_cv_error <- function(Y_hat,y){
  apply( Y_hat , 2, calc_mse, y )
}

#' @rdname kimono name parsing
#' @keywords internal
#' @param string input name
#' @param sep seperator
#' @return vector of mse
parsing_name <- function(string,sep="___"){

  idx <- grepl("___",string) #might be variables like intercept which do not have any prefix

  result <- data.frame('prefix' = as.character(string), "id" = as.character(string), stringsAsFactors = FALSE)
  split_string <- do.call(rbind,strsplit(as.character(string[idx]),'___'))
  result[ idx,1] <- split_string[,1]
  result[ idx,2] <- split_string[,2]

  result
}

#' @rdname kimono  Calculate crossvalidation to estimate lambda1.se & identify potential underpowered features
#' @keywords internal
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param model - input features with prior names attached to features
#' @param intercept boolean
#' @param seed_cv int - remove randomness for crossvalidation
#' @param folds_cv int - defining the amount of folds
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

#' @rdname kimono  Calculate crossvalidation to estimate lambda1.se & identify potential underpowered features
#' @keywords internal
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param model string - which model to train. currently only sparse group lasso tested
#' @param intercept - intercept only model
#' @param seed - seed of model
#' @param ... - forward all other parameters
#' @return edge list for a given input y and x
train_kimono_sgl  <- function(y, x, model = "sparse.grp.lasso", intercept = TRUE, seed_cv = 1234, ...){

  y <- data.table(scale(y))
  x <- data.table(scale(x))

  # estimate best lambda and identify underpowered features
  cv_result   <- calc_cv_sgl(y, x , seed_cv = seed_cv)

  if(length(cv_result)==0) return(c()) #exit function here if all features got excluded

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

  tibble("predictor" = prefix_covariates$id,
         "value" = beta,
         "r_squared" = r_squared,
         "mse" = mse,
         "predictor_layer" = prefix_covariates$prefix
  )
}


#' @rdname kimono Stability selection and summary of multiple runs
#' @keywords internal
#' @param y data.table - feature to predict
#' @param X data.table - input features with prior names attached to features
#' @param nseeds - specifies how many iterations to run
#' @return edge list for a given input y and x, and statistics on multiple runs
stability_selection <- function(y,x,imputed_data,input_data,ADW_calculation, nseeds){

  #Initialize the Seeds and empty Matrix and Vectors for mse and rsquared of length 100
  seeds <- 1:nseeds


  df_values  <- matrix(ncol = nseeds, nrow = ncol(x) + 1 , 0) #+1 accounting for intercept
  df_mse <- rep(NA,nseeds)
  df_r_squared <- rep(NA,nseeds)

  # Groupsparse does not always return all features - therefore we need a vector with all the names
  # and all the relations to calculate the stability selection

  tmp_names <- parsing_name(colnames(x))
  idx <- data.frame('id' = 1:(ncol(x) + 1),
                    'merged' = c('(Intercept)___(Intercept)',colnames(x)),
                    'names' = c('(Intercept)',tmp_names[,2]),
                    'predictor_layer' = c('(Intercept)',tmp_names[,1])
  )

  # Iterate over each of the seeds and set the seed
  for( i in 1:nseeds){

    # Based on the Seed calculate the Cross Validation to determine the best lambda
    # and calculate the best final fit
    fit <- train_kimono_sgl(y,x, seed_cv = seeds[i])
    #fit <- salasso(y, x,imputed_data,input_data,ADW_calculation,set_seed=seeds[i])

    # For some seeds groupsparse does not return a model these have to be skipped
    if(is.null(fit)){
      df_mse[i] <- NA
      df_r_squared[i] <- NA
      df_values[,i] <- NA

      next
    }

    fit <- cbind(fit, 'merged' = paste(fit$predictor_layer,fit$predictor,sep = '___'))
    fit <- merge(idx, fit, by ='merged', all.x = TRUE )

    df_mse[i] <-  fit$mse[1]
    df_r_squared[i] <-  fit$r_squared[1]
    df_values[,i] <-  fit$value
    #reset fit
    fit <- NULL
  }

  #calc

  #return c() if all NA
  if(sum(rowSums(!is.na(df_values))) == 0) return(c())

  #return c() for intercept only model
  if(sum(rowMeans(df_values[idx$names != "(Intercept)",], na.rm = T) != 0) == 0 ) return(c())

  # Dataframe with columns   predictor
  # The value is the frequency of a feature being included in the different seed models
  # Overall R-Squared and MSE are the averages of all R-Squared and MSEs of the different seed models
  # Selected R-Squared and MSE are the averages of the seed models where at least one feature was selected

  tibble(
    'predictor' =  idx$names,
    'mean_value' = rowMeans(df_values, na.rm = T),
    'sd_value' = apply(df_values,1,function(x){sd(x,na.rm = T)}),
    'mean_rsq' = mean(df_r_squared, na.rm = T),
    'sd_rsq' = sd(df_r_squared, na.rm = T),
    'mean_mse' = mean(df_mse, na.rm = T),
    'sd_mse' = sd(df_mse, na.rm = T),
    'sel_freq' = rowSums(df_values != 0 & !is.na(df_values)) / nseeds,
    'predictor_layer' =  idx$predictor_layer
  )
}
