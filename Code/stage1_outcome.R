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
    select(-mu)
    
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

#' Generate Y<stage.suffix> from Normal distribution as a function of covariates
#'
#' @param study.data dataframe; each row corresponds to a ppt
#' @param stage.suffix 
#' @param coefs named vector of numerics; names should correspond to
#' covariates and values should correspond to coefficients
#' @param sigma standard deviation for normal draws
#' @param include.trt.indicators logical indicating whether the treatment indicators 
#' should be included in the data set that is returned. 
#' @param include.mu logical indicating whether the expectation should be included in the data frame as well.
#' Defaults to FALSE because it can cause name collisions with the out of sample data 
#'
#' @return dataframe with Y<suffix> appended. Depending on arguments,
#' may also include A<suffix>_trt inidcator columns or Mu<suffix> expected outcome column
GenNormalOutcomeByFormula <- function(df,
                        stage.suffix,
                        treatment.arm.df,
                        arm.var.name,
                        true.model.formula, 
                        true.param.vec,
                        sigma.noise,
                        include.trt.indicators = TRUE,
                        include.mu = FALSE){
  
  study_data_w_inds <- GenTreatmentIndicators(study.data = df,
                                       treatment.arm.map = treatment.arm.df,
                                       arm.var.name.set = arm.var.name)
  
  full_model_X <- model.matrix(true.model.formula, data=study_data_w_inds)
  
  unused_columns <- setdiff(colnames(full_model_X), names(true.param.vec))
  unused_params <- setdiff(names(true.param.vec), colnames(full_model_X))
  
  # Make sure the order of the columns matches the parameters
  model_X <- full_model_X[, names(true.param.vec)]
  
  # Check for unused columns or parameters and give a warning if there are any
  if (any(c(is_empty(unused_columns), is_empty(unused_params)) == FALSE)) {
    unused_column_string <- paste(unused_columns, collapse = ", ")
    unused_param_string <- paste(unused_params, collapse = ", ")
    
    mismatch_warning_message <- paste0("Column names defined by model formula do not match the parameter names. \n",
                                       "Columns not in parameter vector: ", unused_column_string, "\n",
                                       "Parameters not in covariate matrix: ", unused_param_string)
    warning(mismatch_warning_message)
  }
  
  # Check that the design matrix and the parameter names match
  # Because we reorder the columns this should never happen
  if(all(colnames(model_X) == names(true.param.vec)) == FALSE){
    stop("Covariate matrix columns do not match parameter vector.")
  }
  
  mu_var <- sym(paste0("Mu", stage.suffix))
  outcome_var <- sym(paste0("Y", stage.suffix))
  
  # Calculate the observed outcome (Y<suffix>) and with no noise (Mu<suffix>)
  data_w_outcomes <- study_data_w_inds %>% 
    ungroup %>% 
    mutate(!!mu_var := c(model_X %*% true.param.vec),
           !!outcome_var := !!mu_var + rnorm(n = n(), mean = 0, sd = sigma.noise))
  
  if (include.mu == FALSE) data_w_outcomes <- data_w_outcomes %>% select(-!!mu_var)
  
  if (include.trt.indicators == FALSE) data_w_outcomes <- data_w_outcomes %>% select(-starts_with(paste0("A", stage.suffix)))
  
  return(data_w_outcomes)
}