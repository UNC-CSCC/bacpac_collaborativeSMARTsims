# Functions for allocating first line treatments
# Author: Nikki Freeman
# Last modified on 17 November 2020


#' Randomly assign the first line treatments; no blocking
#'
#' @param df dataframe; each row contains a ppt
#' @param firstLineTreatments character vector of the first line treatments
#'
#' @return dataframe with the first line treatment assignments appended
allocationFn_stage1_v1 <- function(df, firstLineTreatments){
  out <- df %>% rowwise() %>%
    mutate(A1 = sample(x = firstLineTreatments,
                       size = 1)) %>%
    ungroup()
  return(out)
}