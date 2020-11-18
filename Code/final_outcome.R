# Functions for generating Y
# Author: Nikki Freeman
# Last modified on 18 November 2020

#' Add Y1 and Y2 together to make Y
#'
#' @param df dataframe
#'
#' @return dataframe with Y appended where Y = Y1 + Y2
addY1AndY2 <- function(df){
  out <- df %>% mutate(Y = Y1 + Y2)
  
  return(out)
}

#' Set Y to Y2
#'
#' @param df dataframe
#'
#' @return dataframe with a new variable Y where Y = Y2
setYtoY2 <- function(df){
  out <- df %>% mutate(Y = Y2)
  
  return(out)
}

# YforOneDominantRegime <- function(df, stage1TailoringVars, stage1TailoringValues, stage1decision,
#                                   stage2TailoringVars, stage2TailoringValues, stage2decision,
#                                   outcomeDistribution, optRegimeValue, nonOptRegimeValues){
#   nStage1TailoringVars = length(stage1TailoringVars)
#   for(i in 1:nStage1TailoringVars){
#     name = paste("consistent", 1, sep = "_")
#     df <- df %>% mutate(!!name := if_else(!!sym(stage1TailoringVars[i]) == stage1TailoringValues[i], 1, 0))
#   }
#   
#   return(df)
# }


# Optimal regime: If X = 1, then give A1 = "2". Then if response is medium or bad, switch to A2 = "3", if response is good, stay on "2"
# If X= 0, then give A1 = "4". Then if response is medium or bad, switch to "2"; if response good, stay on "4"
YforOneDominantRegime_v1 <- function(df, delta){
  out <- df %>% mutate(consistent = if_else(X_1 == 1 & respStatus == "good" & A1 == "2", 1, 0)) %>%
    mutate(consistent = if_else(X_1 == 1 & A1 == "2" & respStatus != "good" & A2 == "3", 1, consistent)) %>%
    mutate(consistent = if_else(X_1 == 0 & respStatus == "good" & A1 == "4", 1, consistent)) %>%
    mutate(consistent = if_else(X_1 == 0 & A1 == "4" & respStatus != "good" & A2 == "2", 1, consistent)) %>%
    rowwise() %>%
    mutate(Y = if_else(consistent == 1, rnorm(n = 1, mean = delta, sd = 1), rnorm(1, 0, 1))) %>%
    ungroup()
  
  return(out)
}
