# Functions for generating Y2
# Author: Nikki Freeman
# Last modified on 17 November 2020

#' Draw Y2 from a normal distribution
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param coefs named vector of coefficients; coefficients should be the
#' terms of the linear predictor for the mean and the values should be
#' the coefficients
#' @param sigma numeric; standard deviation for the draw from the normal
#'
#' @return dataframe with Y2 appended
generateY2Fn_v1 <- function(df, coefs, sigma){
  out <- df %>% mutate(temp = 1, A1_copy = A1) %>%
    pivot_wider(names_from = A1_copy, names_prefix = "A1_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(temp = 1, A2_copy = A2) %>%
    pivot_wider(names_from = A2_copy, names_prefix = "A2_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(mu := !!(linearPredictor2(intercept = 0, 
                                     coefs = coefs))) %>%
    mutate(mu = Y1 + mu) %>%
    rowwise() %>%
    mutate(Y2 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
    ungroup() %>%
    select(-mu)
  
  return(out)
}
