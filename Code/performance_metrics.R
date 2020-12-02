# Performance Metric Functions
# Author: John Sperger
# Last modified on 01 December 2020


################################################################################
### Design idea:
#
################################################################################

#' Assumes there is a final, overall outcome Y to determine the optimal treatment 
#' sequence. 

calcOptTxAndValFromPOs <- function(mu_df){
  mu_df <- mu_df %>% 
    group_by(ptid) %>% 
    mutate(maxMu = max(Mu),
           OptA1 = A1[which.max(Mu)],
           OptA2 = A2[which.max(Mu)],
           OutcomeUnderOpt = PO[which.max(Mu)],
           DeltaMu = Mu - maxMu,
           DeltaPO = PO - OutcomeUnderOpt) %>% 
    ungroup
  
  return(mu_df)
}

#' Assumes there is a final, overall outcome Y to determine the optimal treatment 
#' sequence. 

summarizeOptTxAndValFromPOs <- function(mu_df){
  opt_summary <- mu_df %>% 
    group_by(ptid) %>% 
    summarise(maxMu = max(Mu),
              OptA1 = A1[which.max(Mu)],
              OptA2 = A2[which.max(Mu)],
              OutcomeUnderOpt = PO[which.max(Mu)])
  
  return(opt_summary)
}

#' Predict the optimal treatment sequence for new observations from a list of fitted
#' models. 
#' @param dtr_fit_both_stages a two element list where the first entry is the fitted 
#' first stage model, and the second entry is the fitted second stage model
#' @param newdata data farme to predict the optimal treatments for
#' @return a n by 3 data frame with columns ptid, A1, and A2
predictTreatSequence <- function(dtr_fit_both_stages,
                                 newdata){
  opt_trt_pred <- tibble(ptid = newdata$ptid,
                     A1 = optTx(dtr_fit_both_stages[[1]], newdata = newdata)$optimalTx,
                     A2 = optTx(dtr_fit_both_stages[[2]], newdata = newdata)$optimalTx)
  
  return(opt_trt_pred)
}


getValueDif <- function(predicted_seq_df,
                        po_grid_df){
  
  #TODO: Because the POs have random noise, the left join will create duplicate rows.
  # It doesn't affect Mu. 
  predicted_seq_with_pos <- predicted_seq_df %>% 
    left_join(x = ., y = po_grid_df, by = c("ptid", "A1", "A2")) %>% 
    distinct_at(., .vars = "ptid", .keep_all = TRUE)
  
  
  return(predicted_seq_with_pos)
}

getValueDifSummary <- function(predicted_seq_df,
                               po_grid_df){
  
  predicted_seq_summary <- getValueDif(predicted_seq_df,
                                       po_grid_df) %>% 
    summarise(MeanVal = mean(Mu),
              PercOracle = mean(Mu)/mean(maxMu))
  
  return(predicted_seq_summary)
}