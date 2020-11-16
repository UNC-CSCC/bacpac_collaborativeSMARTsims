# Functions for BACPAC simulations
# Author: Nikki Freeman
# Last modified on 11 November 2020

# Utility functions -----------------------------------------------------------
#' Helper function to make a linear predictor expression
#'
#' @param intercept integer; interger in linear predictor
#' @param coefs named vector of numerics; names should be the variable name in 
#' the dataframe and the value should be the coefficient
#'
#' @return expression; can be used in the rhs of a mutate statement
linearPredictor <- function(intercept, coefs){
  coef_sym <- syms(names(coefs))
  summands <- map2(coef_sym, coefs, ~ expr((!!.x * !!.y)))
  summands <- c(intercept, summands)
  eq <- reduce(summands, ~ expr(!!.x + !!.y))
  
  return(eq)
}

#' Helper function to make a linear predictor expresssion with interactions
#'
#' @param intercept integer; integer in linear predictor
#' @param coefs named vector of numerics; names should be the variable
#' name(s), e.g., X_1 or X_1*A_2, and the values should be the 
#' corresponding coefficients
#'
#' @return expression; can be used in the rhs of a mutate statement
linearPredictor2 <- function(intercept, coefs){
  # separate the interaction terms and main effects terms
  coefNames <- names(coefs)
  coef_interaction <- coefs[str_detect(coefNames, "\\*")]
  
  if(length(coef_interaction) > 0){
    coef_interaction_sym <- coefNames[str_detect(coefNames, "\\*")]
    coef_main <- coefs[str_detect(coefNames, "\\*", negate = TRUE)]
    coef_main_sym <- syms(names(coef_main))
    
    # Convert the interactions into expressions
    coef_interaction_sym <- str_split(coef_interaction_sym, pattern = "\\*")
    coef_interaction_sym <- map(coef_interaction_sym, linearPredictor2Helper)
    
    # Create summands
    summand_interaction <- map2(coef_interaction_sym, 
                                coef_interaction, ~ expr((!!.x * !!.y)))
    summand_main <- map2(coef_main_sym, coef_main, ~ expr((!!.x * !!.y)))
    # Create the summand
    summands <- c(intercept, summand_interaction, summand_main)
    
    # Create the expression
    eq <- reduce(summands, ~ expr(!!.x + !!.y))
    
    return(eq)
  } else{
    return(linearPredictor(intercept, coefs))
  }
  
  
  return(eq)
}

#' Helper function for linearPredictor2
#'
#' @param coef_interaction_sym character vector of coefficient names
#' for the interaction terms
#'
#' @return coefficient expressions for the interaction terms
linearPredictor2Helper <- function(coef_interaction_sym){
  coef_symbols <- syms(coef_interaction_sym)
  coef_expr <- reduce(coef_symbols, ~ expr(!!.x * !!.y))
  return(coef_expr)
}

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

# Module 1: Make patients and covariates --------------------------------------
#' Make study participants
#'
#' @param N integer; total number of study participants
#'
#' @return dataframe with patient identifiers
makePpts <- function(N){
  out <- data.frame(ptid = 1:N)
  
  return(out)
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
  out <- df %>% exec(covariateFn, !!!..., df = .)
  
  return(out)
}


#' Make covariate function that makes independing bernoulli and normal covars
#'
#' @param df dataframe; each row should correpsond to a ppt 
#' @param numBinaryCovars integer; non-neg number of binary covariates to add
#' @param props vector of probabilities of successes for each binary covariate
#' @param numNormalCovars integer; non-neg number of normal covariates to add
#' @param mu vector of means for the normal covariates
#' @param sd vector of sds for the normal covariates
#'
#' @return dataframe; df with the covariates appended
#' 
#' @export
covariateFn_v1 <- function(df, numBinaryCovars, props,
                           numNormalCovars, mu, sd){
  out <- df %>% 
    makeBinaryCovars(numBinaryCovars, props) %>%
    makeNormalCovariates(numNormalCovars, mu, sd)
  
  return(out)
}
  
#' Add binary covariates to the generated patients
#'
#' @param df dataframe; at a minimum contains the same number of rows 
#' as participants
#' @param numCovars integer; number of binary covariates to add to df
#' @param props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#'
#' @return dataframe with the generated binary covariates appended
makeBinaryCovars <- function(df, numCovars, props){
  if(numCovars < 0){
    warning("Number of binary covariates must be non-negative")
  } else if(numCovars == 0){
    return(df)
  } else{
    covarNames <- paste("X", 1:numCovars, sep = "_")
    for(i in 1:numCovars){
      df <- df %>% rowwise() %>%
        mutate(!!sym(covarNames[i]) := rbinom(n = 1, size = 1, prob = props[i])) %>%
        ungroup()
    }
    return(df)
  }
}


#' Add normal covariates to the generated ppts
#'
#' @param df dataframe; at a minimum contains the same number of rows as
#' participants 
#' @param numCovars non-neg integer; number of normal covars to add
#' @param mu vector; means of normal covariates
#' @param sd vector; sds of normal covariates
#'
#' @return dataframe; df with the normal covariates appended
makeNormalCovariates <- function(df, numCovars, mu, sd){
  if(numCovars < 0){
    warning("Number of normal covariates must be non-negative")
  } else if(numCovars == 0){
    return(df)
  } else{
    covarNames <- paste("W", 1:numCovars, sep = "_")
    for(i in 1:numCovars){
      df <- df %>% rowwise() %>%
        mutate(!!sym(covarNames[i]) := rnorm(n = 1, mean = mu[i], sd = sd[i])) %>%
        ungroup()
    }
    return(df)
  }
}

# Module 2a: Allocate first line treatments -----------------------------------

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

# Module 2b: Generate Y1 and responder status ---------------------------------

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


#' Generate Y1 from Normal distribution as a function of covariates
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param intercept integer; intercept in linear predictor
#' @param coefs named vector of numerics; names should correspond to
#' covariates and values should correspond to coefficients
#' @param sigma standard deviation for normal draws
#'
#' @return dataframe with Y1 appended
generateNormalY1 <- function(df, intercept, coefs, sigma){
  # coef_sym <- syms(names(coefs))
  # summands <- map2(coef_sym, coefs, ~ expr((!!.x * !!.y)))
  # summands <- c(intercept, summands)
  # eq <- reduce(summands, ~ expr(!!.x + !!.y))
  # eq <- linearPredictor(intercept, coefs)

  # out <- df %>% 
  #   mutate(mu := !!eq) %>%
  #   rowwise() %>%
  #   mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
  #   ungroup() %>%
  #   select(-mu)
  
  out <- df %>% mutate(temp = 1, A1_copy = A1) %>%
    pivot_wider(names_from = A1_copy, names_prefix = "A1_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(mu := !!(linearPredictor2(intercept = 0, 
                                     coefs = coefs))) %>%
    mutate(mu = mu + intercept) %>%
    rowwise() %>%
    mutate(Y1 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
    ungroup() %>%
    select(-mu) %>%
    select(-starts_with("A1_"))
  
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

# Module 2c: Allocate second line treatments ----------------------------------

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

# Module 2d. Generate Y2 ------------------------------------------------------
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

#' Draw Y2 from a normal distribution
#'
#' @param df dataframe; each row corresponds to a ppt
#' @param coefs named vector of coefficients; coefficients should be the
#' terms of the linear predictor for the mean and the values should be
#' the coefficients
#' @param sigma numeric; standard deviation for the draw from the normal
#'
#' @return dataframe with Y2 appended
generateY2Fn_v1 <- function(df, coefs, sigma){
  out <- df %>% mutate(temp = 1, A1_copy = A1) %>%
    pivot_wider(names_from = A1_copy, names_prefix = "A1_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(temp = 1, A2_copy = A2) %>%
    pivot_wider(names_from = A2_copy, names_prefix = "A2_", values_from = temp) %>%
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
    mutate(mu := !!(linearPredictor2(intercept = 0, 
                                     coefs = coefs))) %>%
    mutate(mu = Y1 + mu) %>%
    rowwise() %>%
    mutate(Y2 = rnorm(n = 1, mean = mu, sd = sigma)) %>%
    ungroup() %>%
    select(-mu)
  
  return(out)
}

# Module 2e. Generate Y -------------------------------------------------------
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

# Analysis functions ----------------------------------------------------------

doQLearning <- function(dat, moMain_stage2, moCont_stage2,
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
  return(fit_stage1)
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
