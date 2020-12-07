# Utility for BACPAC simulations
# Author: Nikki Freeman
# Last modified on 17 November 2020

# Utility functions -----------------------------------------------------------
#' Helper function to make a linear predictor expression
#'
#' @param intercept integer; interger in linear predictor
#' @param coefs named vector of numerics; names should be the variable name in 
#' the dataframe and the value should be the coefficient
#'
#' @return expression; can be used in the rhs of a mutate statement
linearPredictor <- function(intercept, coefs){
  coef_sym <- syms(names(coefs))
  summands <- map2(coef_sym, coefs, ~ expr((!!.x * !!.y)))
  summands <- c(intercept, summands)
  eq <- reduce(summands, ~ expr(!!.x + !!.y))
  
  return(eq)
}

#' Helper function to make a linear predictor expresssion with interactions
#'
#' @param intercept integer; integer in linear predictor
#' @param coefs named vector of numerics; names should be the variable
#' name(s), e.g., X_1 or X_1*A_2, and the values should be the 
#' corresponding coefficients
#'
#' @return expression; can be used in the rhs of a mutate statement
linearPredictor2 <- function(intercept, coefs){
  # separate the interaction terms and main effects terms
  coefNames <- names(coefs)
  coef_interaction <- coefs[str_detect(coefNames, "\\*")]
  
  if(length(coef_interaction) > 0){
    coef_interaction_sym <- coefNames[str_detect(coefNames, "\\*")]
    coef_main <- coefs[str_detect(coefNames, "\\*", negate = TRUE)]
    coef_main_sym <- syms(names(coef_main))
    
    # Convert the interactions into expressions
    coef_interaction_sym <- str_split(coef_interaction_sym, pattern = "\\*")
    coef_interaction_sym <- map(coef_interaction_sym, linearPredictor2Helper)
    
    # Create summands
    summand_interaction <- map2(coef_interaction_sym, 
                                coef_interaction, ~ expr((!!.x * !!.y)))
    summand_main <- map2(coef_main_sym, coef_main, ~ expr((!!.x * !!.y)))
    # Create the summand
    summands <- c(intercept, summand_interaction, summand_main)
    
    # Create the expression
    eq <- reduce(summands, ~ expr(!!.x + !!.y))
    
    return(eq)
  } else{
    return(linearPredictor(intercept, coefs))
  }
  
  
  return(eq)
}

#' Helper function for linearPredictor2
#'
#' @param coef_interaction_sym character vector of coefficient names
#' for the interaction terms
#'
#' @return coefficient expressions for the interaction terms
linearPredictor2Helper <- function(coef_interaction_sym){
  coef_symbols <- syms(coef_interaction_sym)
  coef_expr <- reduce(coef_symbols, ~ expr(!!.x * !!.y))
  return(coef_expr)
}


#' Generate all possible treatment sequences given vectors of firstline treatment options,
#' second-line treatment options, augmentation treatments (optionally), and a dataframe
#' of impermissible treatment options (optionally). These are all treatment 
#' sequences which could be realized, not the set of feasible treatments for a specific
#' individual. The feasible set of second-line treatments will depend on the observed
#' outcome after stage one. 
#'  
#'  @param first_line_trts character vector of first-line treatments
#'  @param second_line_trts character vector of first-line treatments
#'  @param augment_trts character vector of augmentation treatments
#'  @param impermissible_trt_pairs_df a two column dataframe where the rows are 
#'  impermissible treatment sequences. The first column specifies the first stage 
#'  treatment and the second column specifies the second stage treatment (which
#'  could either be a second-line treatment or an augmentation).
#'  
#'  @return A two column dataframe where the rows are permissible treatment sequences. 
#'  The first column contains the stage one treatmnet, and the second column contains
#'  the second stage treatment (either second-line treatment or augmentation)
generatePossibleTreatmentSequencesGrid <- function(first_line_trts,
                                          second_line_trts,
                                          augment_trts = NULL,
                                          impermissible_trt_pairs_df = NULL){
  
  if(all(first_line_trts %in% second_line_trts) == FALSE) {
    warning("Not all first-line treatments are available in the second-line. 
    Double check that this is what you want. First-line treatments must be
    intentionally specified as an available second-line treatment, otherwise they are not included.")
  }
  
  second_stage_options <- c(second_line_trts, augment_trts)
  
  treatment_pairs_df <- expand_grid(A1 = first_line_trts,
                                    A2 = second_stage_options)
  
  if (is_null(impermissible_trt_pairs_df) == FALSE){
    stopifnot(is.data.frame(impermissible_trt_pairs_df))
    stopifnot(ncol(impermissible_trt_pairs_df) == 2)
    
    if (all(colnames(impermissible_trt_pairs_df) == colnames(treatment_pairs_df)) == FALSE) {
      warning("Column names of impermissible_trt_pairs_df are not A1 and A2,
              the first column will be assumed to be the first stage treatment,
              and the second column the second stage treatment (including augmentations)")
      
      colnames(impermissible_trt_pairs_df) <- colnames(treatment_pairs_df)
    }
    
    #Anti-join returns all rows from x without a match in y
    treatment_pairs_df <- anti_join(x = treatment_pairs_df,
                                    y = impermissible_trt_pairs_df,
                                    by = c("A1", "A2"))
  }
  
  return(treatment_pairs_df)
}