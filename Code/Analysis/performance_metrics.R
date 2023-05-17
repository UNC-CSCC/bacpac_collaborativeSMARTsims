# Performance Metric Functions
# Primary Maintainer: John Sperger 


################################################################################
### Design idea:
#
################################################################################



#' Assumes there is a final, overall outcome Y to determine the optimal treatment 
#' sequence. 

summarizeOptTxAndValFromPOs <- function(mu_df, tailoring_var_names = NULL){
  opt_summary <- mu_df %>% 
    group_by(ptid) %>% 
    summarise(maxMu = max(Mu),
              OptA1 = A1[which.max(Mu)],
              OptA2 = A2[which.max(Mu)],
              OutcomeUnderOpt = PO[which.max(Mu)],
              .groups = 'drop')
  
  if(is_null(tailoring_var_names) == FALSE){
    opt_summary <- mu_df %>% select(ptid, all_of(tailoring_var_names)) %>% 
      distinct_at(., .vars = "ptid", .keep_all = TRUE) %>% 
      right_join(., y = opt_summary, by = "ptid")
  }
  
  return(opt_summary)
}



#' Predict the optimal treatment sequence for new observations from a list of fitted
#' models. 
#' @param dtr_fit_both_stages a two element list where the first entry is the fitted 
#' first stage model, and the second entry is the fitted second stage model
#' @param newdata data farme to predict the optimal treatments for
#' @return a n by 3 data frame with columns ptid, A1, and A2
predictTreatSequence <- function(dtr_fit_both_stages,
                                 newdata,
                                 ...){

  opt_trt_pred <- newdata %>% select(-A1, -A2) %>% 
    mutate(A1 = DynTxRegime::optTx(dtr_fit_both_stages[[1]], newdata = .)$optimalTx)
  
  if ("Y1" %in% colnames(newdata) == FALSE & "Mu1" %in% colnames(newdata) == TRUE) {
    opt_trt_pred <- opt_trt_pred %>% mutate(Y1 = Mu1)
  }
  
  
  opt_trt_pred$A2 <- DynTxRegime::optTx(dtr_fit_both_stages[[2]], 
                                        newdata = as.data.frame(opt_trt_pred))$optimalTx
                     
  
  return(opt_trt_pred %>% select(ptid, A1, A2))
}


getValueDif <- function(predicted_seq_df,
                        po_grid_df,
                        ...){
  
  #TODO: Because the POs have random noise, the left join will create duplicate rows.
  # It doesn't affect Mu. 
  predicted_seq_with_pos <- predicted_seq_df %>% 
    left_join(x = ., y = po_grid_df, by = c("ptid", "A1", "A2")) %>% 
    distinct_at(., .vars = "ptid", .keep_all = TRUE)
  
  
  return(predicted_seq_with_pos)
}

getValueDifSummary <- function(predicted_seq_df,
                               po_grid_df,
                               ...){
  
  predicted_seq_summary <- getValueDif(predicted_seq_df,
                                       po_grid_df) %>% 
    summarise(MeanVal = mean(Mu),
              OracleVal = mean(maxMu),
              PercOracle = mean(Mu)/mean(maxMu),
              .groups = 'drop')
  
  return(predicted_seq_summary)
}
