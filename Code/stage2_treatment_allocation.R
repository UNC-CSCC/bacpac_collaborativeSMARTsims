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


#' Helper function to balance 2nd treatment allocation for those who respond bad
#'
#' @param A1 string; first line treatment assigned
#' @param n integer; number of people assigned to A1
#' @param data dataframe; filtered to those that got A1
#' @param secondLineTreatments character vector; possible second line treatments
#'
#' @return dataframe; data with the second treatment assignments appended
#' @export
#'
#' @examples
allocation2helper1 <- function(A1, n, data, secondLineTreatments){
  feasibleSecondLine <- secondLineTreatments[secondLineTreatments != A1]
  out <- data %>% mutate(A2 = rep_len(gtools::permute(feasibleSecondLine), length.out = n))
  return(out)
}

#' Helper function to balance 2nd treatment allocation for those who respond med and had A1 == SOC
#'
#' @param n integer; number of people assigned to A1 = SOC
#' @param data dataframe; people with A1 = SOC and response = medium
#' @param secondLineTreatments character vector; second line treatments
#' @param augmentationTreatments character vector; augmentation treatments
#' @param standardOfCareTreatment character; standard of care treatment
#'
#' @return dataframe; data wtih A2 appended
#' @export
#'
#' @examples
allocation2helper2 <- function(n, data, secondLineTreatments, augmentationTreatments, standardOfCareTreatment){
  feasibleSecondLine <- c(secondLineTreatments, augmentationTreatments, standardOfCareTreatment)
  out <- data %>% mutate(A2 = rep_len(gtools::permute(feasibleSecondLine), length.out = n))
  return(out)
}

#' Helperfunction to balance 2nd treatment allocation for those who respond med and A1 != SOC
#'
#' @param data dataframe; people with A1 = a1 and response = medium
#' @param n integer; number of people with A1 = a1 and response = medium
#' @param secondLineTreatments character vector; second line treatments
#' @param augmentationTreatments character vector; augmentation treatments
#'
#' @return dataframe; data with A2 appended
#' @export
#'
#' @examples
allocation2helper3 <- function(data, n, secondLineTreatments, augmentationTreatments){
  feasibleSecondLine <- c(secondLineTreatments, augmentationTreatments)
  out <- data %>% mutate(A2 = rep_len(gtools::permute(feasibleSecondLine), length.out = n))
  return(out)
}

#' Allocate stage 2 treatments (equivalent to blocking)
#'
#' @param df dataframe; contain the ppts, first line treatments, and responder status
#' @param firstLineTreatments character vector; firstline treatments
#' @param secondLineTreatments character vector; secondline treatments
#' @param augmentationTreatments character vector; augmentation treatments
#' @param standardOfCareTreatment character; standard of care treatment
#'
#' @return dataframe with A2 appended; Given the feasible treatments after accounting
#' for A1 and responder status, A2 is allocationed evenly across feasible treatments.
#' @export
#'
#' @examples
allocationFn_stage2_v2 <- function(df, firstLineTreatments,
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
    group_by(A1) %>%
    add_count() %>%
    group_by(A1, n) %>%
    nest() %>%
    mutate(data = map(data, allocation2helper1, A1 = A1, n = n, 
                      secondLineTreatments = secondLineTreatments)) %>%
    unnest(data) %>%
    ungroup() %>%
    select(names(good))
    
  medium_soc <- df %>% filter(respStatus == "medium" & A1 == standardOfCareTreatment) %>%
    add_count() %>%
    group_by(A1, n) %>%
    nest() %>%
    mutate(data = map(data, allocation2helper2, n = n, secondLineTreatments = secondLineTreatments,
                      augmentationTreatments = augmentationTreatments, 
                      standardOfCareTreatment = standardOfCareTreatment)) %>%
    unnest(data) %>%
    ungroup() %>%
    select(names(good))
  
  medium <- df %>% filter(respStatus == "medium" & A1 != standardOfCareTreatment) %>%
    group_by(A1) %>%
    add_count() %>%
    group_by(A1, n) %>%
    nest() %>%
    mutate(data = map(data, allocation2helper3, n = n, secondLineTreatments = secondLineTreatments, 
                      augmentationTreatments = augmentationTreatments)) %>%
    unnest(data) %>%
    ungroup() %>%
    select(names(good))
  
  out <- bind_rows(good, bad, medium_soc, medium) %>%
    arrange(ptid)
  
  return(out)
}