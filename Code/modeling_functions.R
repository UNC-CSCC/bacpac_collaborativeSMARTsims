#' Define the feasible set of A2 by A1 but not respStatus for use with DynTxRegime.
#' This makes the assumption that all feasible paths are realized. 
#' It 
#' @param study.data Data frame containing data from the simulated study
#' 
#' @return 
defineFeasibleSet <- function(study.data){
  feasible_A2_by_A1 <- study.data %>% group_by(A1) %>% 
    summarise(Feasible = list(unique(A2)), .groups = "drop_last") %>% 
    mutate(SubgroupName = paste0("subA1_", A1))
  
  subset_list <- map2(.x = feasible_A2_by_A1$SubgroupName,
                      .y = feasible_A2_by_A1$Feasible,
                      ~list(.x, .y))
  
  tx_opts <- left_join(study.data, feasible_A2_by_A1, by = "A1") %>% 
    pull(SubgroupName)
  
  return(list("subsets" = subset_list,
              "txOpts" = tx_opts))
}

#' Define the feasible set of A2 by A1 but not respStatus for use with DynTxRegime.
#' This makes the assumption that all feasible paths are realized. 
#' It 
#' @param study.data Data frame containing data from the simulated study
#' 
#' @return 
fSetFun <- function(A1) {
  if (A1  == 0) {
    subset <- list('subSOC',c(0, 1, 2))
  } else {
    subset <- list('subActive',c(1, 2, 3) )
  }
  return(subset)
}

#' Define the feasible set of A2 by A1 but not respStatus for use with DynTxRegime.
#' This makes the assumption that all feasible paths are realized. 
#' It 
#' @param study.data Data frame containing data from the simulated study
#' 
#' @return 
fSetFunScenario2 <- function(A1) {
  if (A1  == 0) {
    subset <- list('subSOC', c(0, 1, 2, 3))
  } else {
    subset <- list('subActive',c(1, 2, 3) )
  }
  return(subset)
}


#' Define the feasible set of A2 by A1 but not respStatus for use with DynTxRegime.
#' This makes the assumption that all feasible paths are realized. 
#' It 
#' @param study.data Data frame containing data from the simulated study
#' 
#' @return 
fSetFunScenario3 <- function(A1) {
  if (A1  == 0) {
    subset <- list('subSOC', c(0, 1, 2, 3))
  } else {
    subset <- list('subActive',c(1, 2, 3, 4) )
  }
  return(subset)
}
