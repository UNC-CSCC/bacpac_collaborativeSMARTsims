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

#' Randomly assign the first line treatments; equivelent to blocked randomization
#'
#' @param df dataframe; each row contains a ppt
#' @param firstLineTreatments character vector of the first line treatments
#'
#' @return dataframe with the first line treatment assignemnts appended;
#' This function allocates the treatments evenly. There is no need to randomize
#' within the blocks for the simulation; as written this is approximately 
#' the allocation one would get with block randomization
allocationFn_stage1_v2 <- function(df, firstLineTreatments){

  out <- df %>% 
    mutate(A1 = rep_len(firstLineTreatments, length.out = nrow(df)))
  
  return(out)
}


