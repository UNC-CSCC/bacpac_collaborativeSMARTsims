# Functions for Potential Outcomes 
# Author: John Sperger
# Last modified on 18 November 2020

#' Calculate the true mean and potential outcomes for all outcome combinations
#' @param indf dataframe
#' @param widen Logical variable indicating whether a wide dataframe or a long dataframe should be returned. If 
#' wide, the expectations and potential outcomes have names Mu{stage#}_{outcome#}
#' and PO{stage#}_{outcome#}
#' @param trt_options_stage1 vector of treatment options for stage 1
#' @param trt_options_stage2 vector of treatment options for stage 2
#' @param outcome_function_stage1 character string with the name of the function to use for generating Y1
#' @param outcome_function_stage2 character string with the name of the function to use for generating Y2
#' @param outcome_function_final character string with the name of the function to use for generating Y
#' @param assign_responder_fn character string with the name of the function to use for generating repsponder status
#' @param outcome_args_stage1 named list with arguments passed to \code{outcome_function_stage1}
#' @param outcome_args_stage2 named list with arguments passed to \code{outcome_function_stage2}
#' @param outcome_args_final named list with arguments passed to \code{outcome_function_final}
#' @param assign_responder_args named list with arguments passed to \code{assign_responder_fn}
#' @return dataframe with expectations and potential outcomes
calcTrueMeansAndPOs <- function(indf,
                                widen = FALSE,
                                trt_options_stage1,
                                trt_options_stage2,
                                outcome_function_stage1,
                                outcome_function_stage2,
                                outcome_function_final, 
                                outcome_args_stage1,
                                outcome_args_stage2,
                                outcome_args_final,
                                assign_responder_fn,
                                assign_responder_args){
  
  potential_outcomes_long <- .calcTrueMeansAndPOsSingleStage(indf = indf,
                                  stage = 1,
                                  trt_options = trt_options_stage1,
                                  outcome_function = outcome_function_stage1,
                                  outcome_args = outcome_args_stage1) %>% 
    # Temporary workaround because pivot_wider in generateY2Fn_v1 expects unique rows
    mutate(HackyNoiseFix = rnorm(n = length(ptid), mean = 0, sd = 10)) %>% 
    .calcTrueMeansAndPOsSingleStage(indf = .,
                                    stage = 2,
                                    trt_options = trt_options_stage2,
                                    outcome_function = outcome_function_stage2,
                                    outcome_args = outcome_args_stage2) %>% 
    select(-HackyNoiseFix) %>% 
    exec(assign_responder_fn, !!!assign_responder_args, df = .) %>% 
    exec(outcome_function_final, !!!outcome_args_final, df = .) 
    
  if (widen == TRUE){
    indf_with_pos <- .calcTrueMeansAndPOsPivotWide(potential_outcomes_long) %>% 
      full_join(indf, y = ., by = "ptid")
  }
  
  if (widen == FALSE){
    indf_with_pos <- potential_outcomes_long
  }
  
  
  return(indf_with_pos)
}


#' Calculate the true mean for all outcome combinations
#' If the input dataframe contains observed outcomes, the potential outcome will be made to 
#' match the observed outcome. 
#' @param indf dataframe
#' @param stage The stage to calculate the outcome for, either 1 or 2
#' @param trt_options vector of treatment options
#' @param outcome_function character string with the name of the function to use for generating the outcome
#' @param outcome_args named list with arguments passed to \code{outcome_function}
#' @return long dataframe with expectations and potential outcomes potential outcomes
.calcTrueMeansAndPOsSingleStage <- function(indf,
                                            stage,
                                trt_options,
                                outcome_function, 
                                outcome_args){
  
  stopifnot(all(c("coefs", "sigma") %in% names(outcome_args)))
  stopifnot(stage %in% c(1, 2))
  
  
  trt_var_name <- paste0("A", stage)
  observed_outcome_var_name <- paste0("Y", stage)
  expectation_var_name <- paste0("Mu", stage)
  po_var_name <- paste0("PO", stage)
  
  # The oracle knows E[Y | X, A] but not the observed Y | X, A
  # So we'll first calculate mu using no noise, then add the noise back to 
  # get the potential outcomes, to make life easier later
  
  y_sigma <- outcome_args$sigma
  outcome_args$sigma <- 0
  
  po_df <- map(trt_options, 
                            function(in_trt)
                              indf %>%  mutate(!!trt_var_name := in_trt)) %>% 
    bind_rows %>% 
    exec(outcome_function, !!!outcome_args, df = .) %>% 
    rename(!!expectation_var_name := !!sym(observed_outcome_var_name)) %>% 
    mutate(!!po_var_name := !!sym(expectation_var_name) + rnorm(n = 1, mean = 0, sd = y_sigma))
  
  if (observed_outcome_var_name %in% colnames(indf) == FALSE){
    po_df <- po_df %>%  mutate(!!observed_outcome_var_name := !!sym(po_var_name))
    return(po_df)
  }
  
  if (observed_outcome_var_name %in% colnames(indf)){
    # In case you want the potential outcomes/expectations for the study sample,
    # Need to make the potential outcome match the observed outcome for the treatment received
    obs_trt_var_name <- paste0(trt_var_name, "Obs")
    
    indf_plus_pos <- indf %>% 
      select(ptid, all_of(c(trt_var_name, observed_outcome_var_name))) %>% 
      rename(!!obs_trt_var_name := !!trt_var_name) %>% 
      full_join(., po_df, by = "ptid")
      
    for(cur_row in 1:nrow(indf_plus_pos)){
      if(indf_plus_pos[cur_row, obs_trt_var_name] == indf_plus_pos[cur_row, trt_var_name]){
        indf_plus_pos[cur_row, po_var_name] <- indf_plus_pos[cur_row, observed_outcome_var_name] 
      }
    }
    indf_plus_pos <- indf_plus_pos %>% 
      select(-all_of(obs_trt_var_name)) %>% 
      mutate(!!observed_outcome_var_name := !!sym(po_var_name))
    
    return(indf_plus_pos)
  } 
  
  stop("I should be unreachable")
}


.calcTrueMeansAndPOsPivotWide <- function(df_to_pivot){
  
  stage_one_pos_wide_df <- df_to_pivot %>% 
    select(ptid, A1, Mu1, PO1) %>% 
    distinct_at(., .vars = c("ptid", "A1"), .keep_all = TRUE) %>% 
    pivot_wider(names_from = "A1",
                         names_sep = "_",
                         values_from = c("Mu1", "PO1")) 
  
  stage_two_pos_wide_df <- df_to_pivot %>% 
    select(ptid, A1, A2, Mu2, PO2) %>% 
    #distinct_at(., .vars = c("ptid", "A2"), .keep_all = TRUE) %>% 
    pivot_wider(names_from = c("A1", "A2"),
                names_sep = "_",
                values_from = c("Mu2", "PO2"))
  
  final_pos_wide_df <- df_to_pivot %>% 
    select(ptid, A1, A2, Y) %>% 
    pivot_wider(names_from = c("A1", "A2"),
                names_prefix = "Y_",
                values_from = "Y")
  
  pos_wide_df <- full_join(stage_one_pos_wide_df,
                           stage_two_pos_wide_df,
                           by = "ptid") %>% 
    full_join(., final_pos_wide_df, by = "ptid")
  
  return(pos_wide_df)
}