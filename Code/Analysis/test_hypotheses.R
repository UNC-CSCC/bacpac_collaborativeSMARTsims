# Functions related to Inference and Hypothesis Testing
# Primary Maintainer: John Sperger 

################################################################################
### Hypothesis Tests
#
################################################################################

#' 
#' This function uses \code{lm} and 
#' Wald tests to calculate unadjusted p-values.
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
CalcWaldUnadjustedPvals <- function(study.data,
                                          resp.formula,
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

#' Wrapper Around the \code{TukeyHSD} function
#' 
#' 
#' @return vector of *adjusted* p-values for pairwise comparisons

CalcTukeyHSDPvals <- function(study.data,
                         resp.formula,
                         #which.terms,
                         ...){
  
  model_terms <- attr(terms(resp.formula), "term.labels")

  n_model_terms <- length(model_terms)
  stopifnot("Designed to only handle a single model term" = n_model_terms == 1)
  
  # Including second-line treatments gets into some weirdness because A2 depends on Y1 and A1
  # time-varying covariate weirdness, avoiding for now
  
  #stopifnot("Designed to only handle a single model term" = length(which.terms) == 1)
  
  tukey_results <- TukeyHSD(aov(resp.formula, data = study.data))
  
  # TukeyHSD returns a list
  tukey_pvals <- tukey_results[[1]][,4]
  
  return(tukey_pvals)
  
}


#' Define Augment and Switch indicators
#' Use the Wald test to test for equality of the two
TestAugmentVsSwitch <- function(study.data, conditional = FALSE){
  study.data <- study.data %>% 
    mutate(Augment = (A2_0 + A2_1 + A2_2 + A2_3) > 1,
           Switch = (A2 != A1) & Augment == FALSE)
  
  if(conditional == FALSE){
    augment_vs_switch_pvals <- CalcWaldUnadjustedPvals(study.data = study.data %>% filter(A1 != 0),
                                                      resp.formula = formula(Y2 ~ Augment + Switch),
                                                      contrast.mat = matrix(c(0, 1, -1), nrow = 1))
  }
  if(conditional == TRUE){
    switchable_stg1_trts <- 1:3
    augment_vs_switch_pvals <- vector(length = length(switchable_stg1_trts))
    
    for(i in 1:length(switchable_stg1_trts)){
      augment_vs_switch_pvals[i] <- CalcWaldUnadjustedPvals(study.data = study.data %>% filter(A1 == switchable_stg1_trts[i]),
                                                            resp.formula = formula(Y2 ~ Augment + Switch),
                                                            contrast.mat = matrix(c(0, 1, -1), nrow = 1))
    }
    names(augment_vs_switch_pvals) <- paste0("A1_", 1:3)
  }
  
  return(augment_vs_switch_pvals)
}

TestBestNonresponse <- function(study.data, nonresp.cats = c("bad", "medium")){
  nonresp_data <- study.data %>% filter(respStatus %in% nonresp.cats, A2 %in% 0:3) %>% 
    select(A1, A2, Y2)
}
################################################################################
### Adjust P-values
#
################################################################################


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