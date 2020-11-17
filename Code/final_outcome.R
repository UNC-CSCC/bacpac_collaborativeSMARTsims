# Functions for generating Y
# Author: Nikki Freeman
# Last modified on 17 November 2020

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