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

#' Add binary covariates to the generated patients
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
