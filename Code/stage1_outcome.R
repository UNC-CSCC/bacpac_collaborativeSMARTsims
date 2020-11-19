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
  
  # Get the unique treatments
  treatments <- unique(df$A1)
  treatments <- paste("A1", treatments, sep = "_")
  
  # Create placeholders for treatment dummies and add to df
  treatment_df <- data.frame(matrix(0, ncol = length(treatments), nrow = nrow(df)))
  names(treatment_df) <- treatments
  df <- bind_cols(df, treatment_df)
  
  out <- df %>%
    group_by(A1) %>%
    nest() %>%
    mutate(data = map(data, makeDummyHelper_A1, A1 = A1)) %>%
    unnest(data) %>%
    ungroup() %>%
    mutate(mu := !!(linearPredictor2(intercept = 0,
                                     coefs = coefs))) %>%
    mutate(mu = mu + intercept) %>%
    rowwise() %>%
    mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
    ungroup() %>%
    select(-mu) %>%
    select(-starts_with("A1_"))
    
  # out <- simDF %>% mutate(temp = 1, A1_copy = A1) %>%
  #   pivot_wider(names_from = A1_copy, names_prefix = "A1_", values_from = temp) %>%
  #   mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
  #   mutate(mu := !!(linearPredictor2(intercept = 0, 
  #                                    coefs = coefs))) %>%
  #   mutate(mu = mu + intercept) %>%
  #   rowwise() %>%
  #   mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
  #   ungroup() %>%
  #   select(-mu) %>%
  #   select(-starts_with("A1_"))
  
  return(out)
}

makeDummyHelper_A1 <- function(data, A1){
  dummy_name <- paste("A1", A1, sep = "_")
  
  data <- data %>% mutate(!!dummy_name := 1)
  
  return(data)
}
