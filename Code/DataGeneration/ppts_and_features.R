# Functions for making participants and their features 
# Author: Nikki Freeman
# Last modified on 22 November 2020

# Make participants -----------------------------------------------------------
#' Make study participants
#'
#' @param N integer; total number of study participants
#'
#' @return dataframe with patient identifiers
makePpts <- function(N){
  out <- data.frame(ptid = 1:N)
  
  return(out)
}

# Make features/covariates ----------------------------------------------------


#' Make covariate function that makes independing bernoulli and normal covars
#'
#' @param df dataframe; each row should correpsond to a ppt 
#' @param numBinaryCovars integer; non-neg number of binary covariates to add
#' @param props vector of probabilities of successes for each binary covariate
#' @param numNormalCovars integer; non-neg number of normal covariates to add
#' @param mu vector of means for the normal covariates
#' @param sd vector of sds for the normal covariates
#'
#' @return dataframe; df with the covariates appended
#' 
#' @export
covariateFn_v1 <- function(df, numBinaryCovars, props,
                           numNormalCovars, mu, sd){
  out <- df %>% 
    makeBinaryCovars(numBinaryCovars, props) %>%
    makeNormalCovariates(numNormalCovars, mu, sd)
  
  return(out)
}


#' Make covariate function that makes independing bernoulli and normal covars
#'
#' @param df dataframe; each row should correpsond to a ppt 
#' @param numBinaryCovars integer; non-neg number of binary covariates to add
#' @param props vector of probabilities of successes for each binary covariate
#' @param numNormalCovars integer; non-neg number of normal covariates to add
#' @param mu vector of means for the normal covariates
#' @param sd vector of sds for the normal covariates
#'
#' @return dataframe; df with the covariates appended
#' 
#' @export
covariateFn_Mvn <- function(df, numBinaryCovars, props,
                            mu.vec, sigma.mat, d.partial.obs, prop.unobserved){
  out <- df %>% 
    makeBinaryCovars(numBinaryCovars, props) %>%
    MakeMVNCovariates(mu.vec, sigma.mat, d.partial.obs, prop.unobserved)
  
  return(out)
}

#' Add binary covariates to the generated patients that take values in {-1, 1}
#'
#' @param df dataframe; at a minimum contains the same number of rows 
#' as participants
#' @param numCovars integer; number of binary covariates to add to df
#' @param props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#'
#' @return dataframe with the generated binary covariates appended
makeBinaryCovars <- function(df, numCovars, props){
  if(numCovars < 0){
    warning("Number of binary covariates must be non-negative")
  } else if(numCovars == 0){
    return(df)
  } else{
    covarNames <- paste("X", 1:numCovars, sep = "_")
    for(i in 1:numCovars){
      df <- df %>% rowwise() %>%
        mutate(!!sym(covarNames[i]) := 2*rbinom(n = 1, size = 1, prob = props[i])-1) %>%
        ungroup()
    }
    return(df)
  }
}

#' Add binary covariates to the generated patients that takes value in {0, 1}
#'
#' @param df dataframe; at a minimum contains the same number of rows 
#' as participants
#' @param numCovars integer; number of binary covariates to add to df
#' @param props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#'
#' @return dataframe with the generated binary covariates appended
makeBinaryCovarsZeroOne <- function(df, numCovars, props){
  if(numCovars < 0){
    warning("Number of binary covariates must be non-negative")
  } else if(numCovars == 0){
    return(df)
  } else{
    covarNames <- paste("X", 1:numCovars, sep = "_")
    for(i in 1:numCovars){
      df <- df %>% rowwise() %>%
        mutate(!!sym(covarNames[i]) := rbinom(n = 1, size = 1, prob = props[i])) %>%
        ungroup()
    }
    return(df)
  }
}


#' Add normal covariates to the generated ppts
#'
#' @param df dataframe; at a minimum contains the same number of rows as
#' participants 
#' @param numCovars non-neg integer; number of normal covars to add
#' @param mu vector; means of normal covariates
#' @param sd vector; sds of normal covariates
#'
#' @return dataframe; df with the normal covariates appended
makeNormalCovariates <- function(df, numCovars, mu, sd){
  if(numCovars < 0){
    warning("Number of normal covariates must be non-negative")
  } else if(numCovars == 0){
    return(df)
  } else{
    covarNames <- paste("W", 1:numCovars, sep = "_")
    for(i in 1:numCovars){
      df <- df %>% rowwise() %>%
        mutate(!!sym(covarNames[i]) := rnorm(n = 1, mean = mu[i], sd = sd[i])) %>%
        ungroup()
    }
    return(df)
  }
}


#' Add multivariate normal covariates to the generated ppts
#' Optionally have some of the covaraites be partially observed. 
#' Creates indicators of whether the covariate was observed for the partially observed covariates
#' @param df dataframe; at a minimum contains the same number of rows as
#' participants 
#' @param mu.vec vector; means of normal covariates
#' @param sd sigma.mat; sds of normal covariates
#'
#' @return dataframe; df with the normal covariates appended. 
#' W_x for continuous covariates
#' U_x for partially observed covariates
#' O_x for indicator of whether the subject had expensive covariates measured
MakeMVNCovariates <- function(df, mu.vec, sigma.mat, d.partial.obs = 0, prop.unobserved){
  if (is.null(mu.vec)) return(df)
  
  cont_covar <- mvtnorm::rmvnorm(n = nrow(df), mean = mu.vec, sigma = sigma.mat)
  
  if (d.partial.obs > 0){
    d.observed <- ncol(cont_covar) - d.partial.obs
    d.true.params <- ncol(cont_covar)
    
    # Want to keep the partially observed variables for data generation
    obs_var_names <- paste0("W_", 1:length(mu.vec))
    partial_var_names <- paste0("U_", 1:d.partial.obs)
    obs_indicator_var_names <- paste0("O_", 1:d.partial.obs)
    
    start_pos_partial_covar <- min(d.true.params, d.true.params - d.partial.obs + 1)
    
    partial_covar <- matrix(cont_covar[, start_pos_partial_covar:d.true.params], ncol = d.partial.obs)
    partial_covar[sample(1:nrow(partial_covar), 
                         size = floor((prop.unobserved)*nrow(partial_covar))), ] <- NA

    obs_indicators <- is.na(partial_covar) == FALSE
      
    cont_covar <- cbind(cont_covar, partial_covar, obs_indicators)
    
    
    colnames(cont_covar) <- c(obs_var_names, partial_var_names, obs_indicator_var_names)
    
    
  } else{ colnames(cont_covar) <- paste0("W_", 1:ncol(cont_covar))}
  
  df <- bind_cols(df, as_tibble(cont_covar))
  
  return(df)
}