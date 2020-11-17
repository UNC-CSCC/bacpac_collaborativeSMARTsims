# Functions assigning responder status
# Author: Nikki Freeman
# Last modified on 17 November 2020

#' Assign responder status by quantile
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param ctsOutcome1 character; name of the cts cut-off variable
#' @param cutoffs vector of quantiles (e.g. c(0.2, 0.8))
#'
#' @return dataframe wtih the responder stsatus appended
assignResponderStatusByQuantile <- function(df, ctsOutcome1, cutoffs){
  out <- df %>%
    group_by(A1) %>%
    mutate(pctRank = percent_rank(!!ensym(ctsOutcome1))) %>%
    mutate(respStatus = if_else(pctRank <= cutoffs[1], "bad", "medium"),
           respStatus = if_else(pctRank >= cutoffs[2], 
                                "good", respStatus))
  return(out)
}