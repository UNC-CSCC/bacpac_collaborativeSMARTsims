# Functions for allocating stage 2 treatments
# Author: Nikki Freeman
# Last modified on 17 November 2020

#' Allocation function for stage 2 that is based on the current schematic (11/11)
#'
#' @param df dataframe; each row corresponds to participant and should include
#' the first line treatments and responder statuses
#' @param firstLineTreatments character vector of first line treatments
#' @param secondLineTreatments character vector of second line treatments (does not 
#' include the standard of care treatment)
#' @param augmentationTreatments character vector of the augmentation treatments
#' @param standardOfCareTreatment character vector with the name of the standard of 
#' care treatment
#'
#' @return dataframe with the second line treatments
allocationFn_stage2_v1 <- function(df, firstLineTreatments,
                                   secondLineTreatments,
                                   augmentationTreatments,
                                   standardOfCareTreatment){
  # Available treatments depend on the response status
  # Consider each response group separately
  # Good is the easiest case, just continue on with A1
  good <- df %>% filter(respStatus == "good") %>%
    mutate(A2 = A1)
  
  # Bad is the second easiest case, switch off of A1 and no augmentation
  bad <- df %>% filter(respStatus == "bad") %>%
    rowwise() %>%
    mutate(A2 = sample(x = secondLineTreatments[secondLineTreatments != A1],
                       size = 1)) %>%
    ungroup()
  
  medium_soc <- df %>% filter(respStatus == "medium" 
                              & A1 == standardOfCareTreatment) %>%
    rowwise() %>%
    mutate(A2 = 
             sample(x = c(secondLineTreatments, 
                          standardOfCareTreatment, 
                          augmentationTreatments),
                    size = 1)) %>%
    ungroup()
  
  medium <- df %>% filter(respStatus == "medium" 
                          & A1 != standardOfCareTreatment) %>%
    rowwise() %>%
    mutate(A2 = sample(x = c(secondLineTreatments, augmentationTreatments),
                       size = 1)) %>%
    ungroup()
  
  out <- bind_rows(good, bad, medium_soc, medium) %>%
    arrange(ptid)
  
  return(out)
}