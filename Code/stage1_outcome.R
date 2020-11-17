# Functions for generating Y1
# Author: Nikki Freeman
# Last modified on 17 November 2020

#' Generate Y1 from Normal distribution as a function of covariates
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param intercept integer; intercept in linear predictor
#' @param coefs named vector of numerics; names should correspond to
#' covariates and values should correspond to coefficients
#' @param sigma standard deviation for normal draws
#'
#' @return dataframe with Y1 appended
generateNormalY1 <- function(df, intercept, coefs, sigma){
  # coef_sym <- syms(names(coefs))
  # summands <- map2(coef_sym, coefs, ~ expr((!!.x * !!.y)))
  # summands <- c(intercept, summands)
  # eq <- reduce(summands, ~ expr(!!.x + !!.y))
  # eq <- linearPredictor(intercept, coefs)
  
  # out <- df %>% 
  #   mutate(mu := !!eq) %>%
  #   rowwise() %>%
  #   mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
  #   ungroup() %>%
  #   select(-mu)
  
  out <- df %>% mutate(temp = 1, A1_copy = A1) %>%
    pivot_wider(names_from = A1_copy, names_prefix = "A1_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(mu := !!(linearPredictor2(intercept = 0, 
                                     coefs = coefs))) %>%
    mutate(mu = mu + intercept) %>%
    rowwise() %>%
    mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
    ungroup() %>%
    select(-mu) %>%
    select(-starts_with("A1_"))
  
  return(out)
}


