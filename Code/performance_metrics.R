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

################################################################################
### Hypothesis Testing (Currently based only on the first stage)
#
################################################################################

#' 
#' This function uses \code{lm}
#' 
#' @param study.data data frame containing the study data
#' @param resp.formula formula object with the linear model to the data
#' @param contrast.mat matrix where each row is a contrast vector for a single df Wald test. Only one of contrast.mat
#' or coefs.names.to.test should be specified
#' @param coefs.names.to.test Character vector where each entry is the name of a parameter to test. This will create the required
#' contrast matrix, but is only capable of testing whether single parameters are equal to zero. 
#' 
#' @return vector of p-values with length = nrow(contrast.mat) of unadjusted Wald p-values
#' 
CalcStage1WaldUnadjustedPvals <- function(study.data,
                          resp.formula = formula(Y1 ~ A1*X_1),
                          contrast.mat = NULL,
                          coefs.to.test = NULL){
  if (is_null(contrast.mat) == FALSE & is_null(coefs.to.test) == FALSE) stop("Only one of contrast.mat and coefs.to.test should be specified at a time")
  if (is_null(contrast.mat) == TRUE & is_null(coefs.to.test) == TRUE) stop("One of contrast.mat and coefs.to.test must be specified")
  
  lm_fit <- lm(formula = resp.formula, data = study.data)
  
  p_params <- length(coef(lm_fit))
  
  if (is_null(coefs.to.test) == FALSE) {
    m_hypotheses <- length(coefs.to.test)
    
    contrast.mat <- matrix(data = 0, nrow = m_hypotheses, ncol = p_params)
    rownames(contrast.mat) <- paste0("Placeholder", 1:m_hypotheses)
    
    for(cur_index in 1:m_hypotheses){
      coef_index <- which(names(coef(lm_fit)) == coefs.to.test[cur_index])
      
      contrast.mat[cur_index, coef_index] <- 1L
      
      rownames(contrast.mat)[cur_index] <- coefs.to.test[cur_index]
    }
  }
  
  Sigma_hat <- vcov(lm_fit)
  
  wald_stats <- vector(length = nrow(contrast.mat))
  wald_pvals <- vector(length = nrow(contrast.mat))
  
  beta_hat <- coef(lm_fit)
  
  for(i in 1:length(wald_stats)){
    current.contrast <- matrix(contrast.mat[i,], nrow=1)
    wald.outer <- (current.contrast %*% beta_hat)
    wald.inner <- solve(current.contrast %*% Sigma_hat %*% t(current.contrast))
    wald_stats[i] <- t(wald.outer) %*% wald.inner %*% wald.outer
    wald_pvals[i] <- pchisq(q = wald_stats[i], df = 1, lower.tail = FALSE)
  }
  
  names(wald_pvals) <- rownames(contrast.mat)
  
  return(wald_pvals)
}

#' Wrapper around p.adjust to adjust p-values for multiple comparisons.  
#' 
#' @param pval.mat An matrix of p-values. The rows should be simulation runs, and 
#' the columns unadjusted p-values for hypotheses.
#' @param adj.method string specifying the adjustment method to use. Defaults to "holm"
#' 
#' @return A matrix of adjusted p=valuues with the same dimensions as pval.mat 
AdjustPvals <- function(pval.mat, adj.method = "holm"){
  adjusted_pval_mat <- t(apply(pval.mat, 1, FUN = p.adjust, method = adj.method))
  colnames(adjusted_pval_mat) <- colnames(pval.mat)
  return(as_tibble(adjusted_pval_mat))
}
