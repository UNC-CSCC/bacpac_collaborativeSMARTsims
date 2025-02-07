# Functions assigning responder status
# Author: Nikki Freeman
# Last modified on 17 November 2020

#' Assign responder status by quantile
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param ctsOutcome1 character; name of the cts cut-off variable
#' @param cutoffs vector of quantiles (e.g. c(0.2, 0.8))
#'
#' @return dataframe wtih the responder status appended
assignResponderStatusByQuantile <- function(df, ctsOutcome1, cutoffs, status.labels = c("bad", "medium", "good")){
  out <- df %>%
    group_by(A1) %>%
    mutate(pctRank = percent_rank(!!ensym(ctsOutcome1))) %>%
    mutate(respStatus = cut(pctRank, 
                            breaks = c(0, cutoffs, 1), 
                            include.lowest = TRUE, 
                            labels = status.labels))
  return(out)
}

#' Assign responder status by quantile from the N(0, 1)
#'
#' @param df dataframe; each row corresponds to a ppts and has Y1
#' @param cutoffs vector of quantiles (e.g. c(0.2, 0.8))
#'
#' @return dataframe with responder status appended
#' @export
#'
#' @examples
assignResponderStatusByNormalQuantileForStdNormalY1 <- function(df, cutoffs,
                                                                status.labels = c("bad", "medium", "good")){
  out <- df %>% 
    group_by(A1) %>%
    mutate(Y1_std = (Y1 - mean(Y1))/sd(Y1)) %>%
    mutate(respStatus =  cut(Y1_std, 
                             breaks = c(-Inf, qnorm(cutoffs), Inf), 
                             labels = status.labels)) %>%
    select(-Y1_std)
  
  return(out)
}

