# Key wrapper functions for simulations
# Author: Nikki Freeman
# Last modified on 17 November 2020

# Functions for making the simulation data ------------------------------------

#' Make one simulation data set
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#'
#' @return data frame with 1 simulation data set
#' @export
make1SimDataSet <- function(metadata, args){
  out <- exec(args$makePpts_fn, !!!args$makePpts_args) %>%
    exec(args$makeCovariates_fn, !!!args$makeCovariates_args, df = .) %>%
    exec(args$allocateStage1Treatments_fn, !!!args$allocateStage1Treatments_args, df = .) %>%
    exec(args$generateY1_fn, !!!args$generateY1_args, df = .) %>%
    exec(args$assignResponderStatus_fn, !!!args$assignResponderStatus_args, df = .) %>%
    exec(args$allocateStage2Treatments_fn, !!!args$allocateStage2Treatments_args, df = .) %>%
    exec(args$generateY2_fn, !!!args$generateY2_args, df = .) %>%
    exec(args$generateY_fn, !!!args$generateY_args, df = .) 
  
  return(out)
}

#' Wrapper function for making simulation data
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#' @param set integer; simulation data set counter
#'
#' @return data frame with 1 simulation data set labeled by dataSet = set
makeSimDataWrapper <- function(metadata, args, set){
  out <- make1SimDataSet(metadata, args) %>%
    mutate(dataSet = set)
  
  return(out)
} 


# Analysis functions ----------------------------------------------------------

doQLearning <- function(dat, moMain_stage2, moCont_stage2,
                        moMain_stage1, moCont_stage1){
  
  fit_stage1 <- doQLearningListReturn(dat = dat,
                        moMain_stage2 = moMain_stage2, 
                        moCont_stage2 = moCont_stage2,
                        moMain_stage1 = moMain_stage1, 
                        moCont_stage1 = moCont_stage1)[[1]]
  
  return(fit_stage1)
}

doQLearningListReturn <- function(dat, moMain_stage2, moCont_stage2,
                        moMain_stage1, moCont_stage1){
  
  dat$A2 <- factor(dat$A2)
  dat <- data.frame(dat)
  
  
  fit_stage2 <- DynTxRegime::qLearn(moMain = moMain_stage2,
                                    moCont = moCont_stage2,
                                    response = dat$Y,
                                    data = dat,
                                    txName = "A2")
  fit_stage1 <- DynTxRegime::qLearn(moMain = moMain_stage1,
                                    moCont = moCont_stage1,
                                    response = fit_stage2,
                                    data = dat,
                                    txName = "A1")
  return(list(stage1 = fit_stage1,
              stage2 = fit_stage2))
}


# Graveyard -------------------------------------------------------------------
if(FALSE){
  #' Make the potential outcomes for stage 1
  #'
  #' @param df dataframe; at a minimum contains one row for each ppt
  #' @param treatment_list list; first element is the number of treatments, 
  #' second element is the name of each treatment; third element is themean of 
  #' each treatment, fourth element is the sd of each treatment; fourth element 
  #' is the randomization prob for each treatment
  #' @param cf_name character; name to give the vector of counterfactuals
  #'
  #' @return dataframe with the counterfactual outcomes for each treatment
  #' where outcomes are drawn from the normal(mu, sigma)
  makeStage1PotentialOutcomes <- function(df, treatment_list, cf_name){
    out <- df %>%
      group_by(ptid) %>%
      nest() %>%
      mutate(!!cf_name := makeCFHelper(data, treatment_list)) %>%
      select(-data)
    
    
    return(out)
    
  }
  
  #' Helper function for making counterfactuals
  #'
  #' @param data df; consequence of nested dataframe
  #' @param treatment_list list; contains the treatment information
  #' including the mean and sd for each treatment
  #'
  #' @return list; named list with the counterfactual outcomes
  makeCFHelper <- function(data, treatment_list){
    out <- map2(treatment_list[["mu_treatments"]], 
                treatment_list[["sd_treatments"]],
                rnorm, n = 1) %>%
      unlist()
    names(out) <- treatment_list[["treatments"]]
    out <- list(out)
    return(out)
  }
  
  
  #' Assign treatments
  #'
  #' @param df dataframe; at a minimum number of rows is equal to the number of participants
  #' @param treatment_list list; contains treatment information for decision point (given 
  #' responder status if applicable)
  #' @param assignment character; name for the treatment assignment, probably A1 or A2
  allocateTreatments <- function(df, treatment_list, assignment){
    out <- df %>%
      rowwise() %>%
      mutate(!!assignment := sample(x = treatment_list[["treatments"]], 
                                    size = 1, 
                                    prob = treatment_list[["probs"]])) 
    return(out)
  }
  
  getSecondStageTreatmentList <- function(A1, respStatus, 
                                          treatment_list,
                                          standardOfCareTreatment,
                                          firstLineAvailTreatmentsList,
                                          # augmentableTreatments, 
                                          augumentationTreatmentsList,
                                          newTreatmentsList){
    if(respStatus == "good"){
      out <- list(n_treatments = 1,
                  treatments = A1,
                  mu_treatments = treatment_list[["mu_treatments"]][treatment_list[["treatments"]] == A1],
                  sd_treatments = treatment_list[["sd_treatments"]][treatment_list[["treatments"]] == A1],
                  probs = 1)
      return(out)
    } else if(respStatus == "bad"){
      # Available treatments are firstline treatments available at the second randomization
      # Additional treatments may also be available at the second randomization
      treatments <- map(list(firstLineAvailTreatmentsList, newTreatmentsList), 
                        pluck, "treatments") %>%
        unlist()
      mu_treatments <- map(list(firstLineAvailTreatmentsList, newTreatmentsList), 
                           pluck, "mu_treatments") %>%
        unlist()
      sd_treatments <- map(list(firstLineAvailTreatmentsList, newTreatmentsList),
                           pluck, "sd_treatments") %>%
        unlist()
      # Secondline treatment cannot be the firstline treatment if doing badly
      mu_treatments <- mu_treatments[treatments != A1]
      sd_treatments <- sd_treatments[treatments != A1]
      treatments <- treatments[treatments != A1]
      
      # Create a list for treatment options
      n_treatments <- length(treatments)
      out <- list(n_treatments = n_treatments,
                  treatments = treatments,
                  mu_treatments = mu_treatments,
                  sd_treatments = sd_treatments,
                  probs = rep(1/n_treatments, n_treatments))
      return(out)
    } else if(respStatus == "medium"){
      # Available treatments are firstline treatments available as secondline treatments
      # Augmentation treatments are available
      # Additional secondline treatments may also be available
      treatments <- map(list(firstLineAvailTreatmentsList, 
                             newTreatmentsList
                             # , 
                             # augmentationTreatmentList
      ), 
      pluck, "treatments") %>%
        unlist() 
      mu_treatments <- map(list(firstLineAvailTreatmentsList, 
                                newTreatmentsList
                                # , 
                                # augmentationTreatmentList
      ), 
      pluck, "mu_treatments") %>%
        unlist()
      sd_treatments <- map(list(firstLineAvailTreatmentsList, 
                                newTreatmentsList
                                # , 
                                # augmentationTreatmentList
      ), 
      pluck, "sd_treatments") %>%
        unlist()
      
      if(A1 == standardOfCareTreatment[["treatments"]]){
        treatments <- c(treatments, pluck(standardOfCareTreatment, "treatments"))
        mu_treatments <- c(mu_treatments, pluck(standardOfCareTreatment, "mu_treatments"))
        sd_treatments <- c(sd_treatments, pluck(standardOfCareTreatment, "sd_treatments"))
      }
      n_treatments <- length(treatments)
      out <- list(n_treatments = n_treatments,
                  treatments = treatments,
                  mu_treatments = mu_treatments,
                  sd_treatments = sd_treatments,
                  probs = rep(1/n_treatments, n_treatments))
      return(out)
    }
  }
  
  
  makeStage2PotentialOutcomes <- function(data){
    out <-  map2(data$treatment_list2[[1]][["mu_treatments"]],
                 data$treatment_list2[[1]][["sd_treatments"]],
                 rnorm, n = 1) %>%
      unlist()
    names(out) <-data$treatment_list2[[1]][["treatments"]]
    return(out)
  }
}

#' Make Covariates
#'
#' @param df dataframe; at a minimum each row should correspond to a ppt
#' @param covariateFn function; function for creating covariates
#' @param ... varargs to pass to covariateFn
#'
#' @return
#' @export
#'
#' @examples
makeCovariates <- function(df, covariateFn, ...){
  out <- df %>% covariateFn(...)
  
  return(out)
}

#' Allocate the stage 1 treatments
#'
#' @param df dataframe; each row should correspond to a ppt
#' @param allocationFn function; allocates the first line treatments
#' @param ... varargs to pass to allocationFn
#'
#' @return dataframe with the allocated treatments appended
allocateStage1Treatments <- function(df, allocationFn, ...){
  out <- df %>% allocationFn(...)
  return(out)
}

#' Generate Y1
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param generateY1Fn function; function that generates Y1
#' @param ... varargs to pass to generateY1Fn
#'
#' @return dataframe with Y1 appended
generateY1 <- function(df, generateY1Fn, ...){
  out <- df %>% generateY1Fn(...)
  
  return(out)
}


#' Assign responder status
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param assignResponderStatusFn function to assign responder status function
#' @param ... varargs to pass to assignResponderStatusFn
#'
#' @return dataframe with the responder status appended
assignResponderStatus <- function(df, assignResponderStatusFn, ...){
  out <- df %>% assignResponderStatusFn(...)
  
  return(out)
}

#' Allocate stage 2 treatments
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param allocationFn_stage2 function to allocate stage 2 treatments
#' @param ... varargs to pass to allocationFn_stage2
#'
#' @return dataframe with the stage 2 treatment allocations
allocateStage2Treatments <- function(df, allocationFn_stage2, ...){
  
  out <- df %>% allocationFn_stage2(...)
  
  return(out)
}

#' Generate Y2
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param generateY2Fn function to generate Y2
#' @param ... varargs to pass to generateY2Fn
#'
#' @return dataframe with Y2 appended
generateY2 <- function(df, generateY2Fn, ...){
  out <- df %>% generateY2Fn(...)
  
  return(out)
}

#' Generate Y
#'
#' @param df data frame
#' @param generateYFn function to generate Y
#' @param ... varargs to pass to generateYFn
#'
#' @return dataframe with Y appended
generateY <- function(df, generateYFn, ...){
  out <- df %>% generateYFn(...)
  
  return(out)
}
