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
                        moMain_stage1, moCont_stage1, fSetFunction = NULL){
  
  fit_stage1 <- doQLearningListReturn(dat = dat,
                        moMain_stage2 = moMain_stage2, 
                        moCont_stage2 = moCont_stage2,
                        moMain_stage1 = moMain_stage1, 
                        moCont_stage1 = moCont_stage1,
                        fSetFunction = fSetFunction)[[1]]
  
  return(fit_stage1)
}

doQLearningListReturn <- function(dat, moMain_stage2, moCont_stage2,
                        moMain_stage1, moCont_stage1, fSetFunction){
  
  dat$A2 <- factor(dat$A2)
  dat$A1 <- factor(dat$A1)
  
  dat <- data.frame(dat)
  
  
  fit_stage2 <- DynTxRegime::qLearn(moMain = moMain_stage2,
                                    moCont = moCont_stage2,
                                    response = dat$Y,
                                    data = dat,
                                    txName = "A2",
                                    fSet = fSetFunction,
                                    verbose = FALSE)
  
  fit_stage1 <- DynTxRegime::qLearn(moMain = moMain_stage1,
                                    moCont = moCont_stage1,
                                    response = fit_stage2,
                                    data = dat,
                                    txName = "A1",
                                    verbose = FALSE)
  return(list(stage1 = fit_stage1,
              stage2 = fit_stage2))
}

#' Make one out-of-sample data set. Currently (as of 12/6/20) this will be a long data
#' frame with a row for each out-of-sample patient and hypothetical treatment sequence
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#'
#' @return long data frame for one out of sample data set 
#' @export 

make1OutofSampleSimDataSet <- function(metadata, args){
 
  out <- exec(args$makePpts_fn, !!!args$makePpts_args) %>%
    exec(args$makeCovariates_fn, !!!args$makeCovariates_args, df = .) %>% 
    exec(args$MakePOs_fn, !!!args$makePOs_args, indf = .)
  

  return(out)
}

#' Create in-sample and out-of-sample performance metrics
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#' @param study.data.list list where each element is a data frame of simulated study data
#' @param oos.data 
#'
#' @return long data frame for one out of sample data set 
#' @export 
analyzeSimulationRunsOneSettingWrapper <- function(metadata, args,
                                                   study.data.list,
                                                   oos.data){
  stop("This function is not finished")
  
  predicted_sequences <- map(study.data.list, 
                             ~exec(args$AnalysisModeling_fn, 
                                   args$AnalysisModeling_args, dat = .)) %>% 
    map(., 
        ~exec(args$PredictSequence_fn, args$PredictSequence_args, newdata = oos.data))

  if (is_empty(args$makeInSampleMetrics_fns) == FALSE){
    # Each entry of args$makeInSampleMetrics_fns is a function for calculating some
    # in-sample performance metric
    
    k_in_sample_metrics <- length(args$makeInSampleMetrics_fns)
    
    in_sample_metric_results <- vector("list", length = k_in_sample_metrics)
    
    for(cur_index in 1:k_in_sample_metrics) {
      in_sample_metric_results[[cur_index]] <- map(study.data.list, 
                                                   ~exec(args$makeInSampleMetrics_fns_list[[cur_index]], 
                                                         !!!makeInSampleMetrics_args_list[[cur_index]],
                                                         indf = .))
    }
  }
  
  if (is_empty(args$makeOOSMetrics_fns) == FALSE){
    # Each entry of args$makeInSampleMetrics_fns is a function for calculating some
    # in-sample performance metric
    
    k_oos_sample_metrics <- length(args$makeOOSMetrics_fns)
    
    oos_metric_results <- vector("list", length = k_oos_sample_metrics)
    
    for(cur_index in 1:k_oos_sample_metrics) {
      oos_metric_results[[cur_index]] <- map(study.data.list, 
                                                   ~exec(args$makeOOSMetrics_fns[[cur_index]], 
                                                         !!!makeOOSMetrics_args[[cur_index]],
                                                         oos.data = oos.data,
                                                         indf = .))
    }
  }
  
  results_list <- c(in_sample_metric_results, oos_metric_results)
  
  return(results_list)
  
}

#' Calculate percentage of oracle value
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#'
#' @return long data frame for one out of sample data set 
#' @export 
analyzeSimulationRunsPercOracleWrapper <- function(metadata, args,
                                                   study.data.list,
                                                   oos.data){
  oos_unique <- oos.data %>% distinct_at(., .vars = "ptid", .keep_all = TRUE)
  
  oos_vals <- map(study.data.list, 
                             ~exec(args$AnalysisModeling_fn, 
                                   !!!args$AnalysisModeling_args, dat = .)) %>% 
    map(., 
        ~exec(args$PredictSequence_fn, 
              dtr_fit_both_stages = .,
              newdata = oos_unique),
              !!!args$PredictSequence_args) %>% 
     map_dfr(., ~getValueDifSummary(., po_grid_df = oos.data)) %>% 
    mutate(N = map_int(study.data.list, nrow))
  
  return(oos_vals)
}

#' Calculate percentage of oracle value
#'
#' @param metadata list of metadata (N, firstline, secondline, augementation, and SOC treatments)
#' @param args list of functions and arguments for each key trial design and simulation aspect
#'
#' @return long data frame for one out of sample data set 
#' @export 
analyzeSimulationRunsPercOracleWrapperParallel <- function(metadata, args,
                                                   study.data.list,
                                                   oos.data){
  
  oos_vals <- furrr::future_map(study.data.list, 
                  ~exec(args$AnalysisModeling_fn, 
                        !!!args$AnalysisModeling_args, dat = .),
                  .options = furrr::furrr_options(seed = TRUE)) %>% 
    furrr::future_map(., 
        ~exec(args$PredictSequence_fn, 
              dtr_fit_both_stages = .,
              newdata = oos.data,
              !!!args$PredictSequence_args),
        .options = furrr::furrr_options(seed = TRUE)) %>% 
    furrr::future_map_dfr(., ~getValueDifSummary(., po_grid_df = oos.data),
                          .options = furrr::furrr_options(seed = TRUE)) %>% 
    mutate(N = map_int(study.data.list, nrow))
  
  return(oos_vals)
}


# Specific Step Wrappers ----------------------------------------------------------

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
