# Functions for allocating first line treatments
# Author: Nikki Freeman
# Last modified on 17 November 2020


################################################################################
## Block randomization
################################################################################

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

#' Stratified Permuted block randomization for stage one 
#' @param df dataframe; each row contains a ppt
#' @param firstLineTreatments character vector of the first line treatments
#' @param strata.var.syms a vector with the variable names to stratify by. The 
#' variable names should be symbols, not strings, for example sym("X_1") instead of "X_1"
#' @return dataframe with the first line treatment appended (column name A1)
StratifiedBlockRandomizationStage1 <- function(df, first.line.trts, 
                                         strata.vars.syms = c(sym("X_1"))) {
  data_w_trts <- df %>% group_by(!!!strata.vars.syms) %>% 
    mutate(A1 = sample(rep(first.line.trts, 
                             each = ceiling(n()/length(first.line.trts))),
                         size = n())) %>% 
    ungroup 
  
  return(data_w_trts)
}

################################################################################
## Minimization
################################################################################


#' 
allocationFn_stage1_Minimization <- function(df, covar.names.to.min, firstLineTreatments, ...){
  # our binary covariates use {-1, 1} coding
  # balanceRandomize() expects factor levels to be 1, 2, ...
  # transforming them to factors and back to numeric achieves this
  covar_matrix <- df %>% 
    select(all_of(covar.names.to.min)) %>% 
    mutate_at(all_of(covar.names.to.min), factor) %>% 
    mutate_at(all_of(covar.names.to.min), as.numeric) %>% 
    as.matrix
  
  allocation <- balanceRandomize(Z = covar_matrix, N = length(firstLineTreatments), ...)

  out <- df %>% mutate(A1 = firstLineTreatments[allocation$assignment])
  
  return(out)
}

#' Applies the minimization algorithm to each cluster independently
#' 
#' @param df
#' @param covar.names.to.min vector of strings with the names of the covariates to include in the minimization
#' @param firstLineTreatments vector of 
AllocateStage1TrtMinCluster <- function(df, covar.names.to.min, firstLineTreatments, ...) {
  stopifnot("Site" %in% colnames(df))
  
  out <-  df %>% 
    group_by(Site) %>% 
    group_modify(., ~allocationFn_stage1_Minimization(.,
                                                      covar.names.to.min = covar.names.to.min, 
                                                      firstLineTreatments = firstLineTreatments)) %>% 
    ungroup
  
  return(out)
}

#' Minimization algorithm based on Pocock and Simon 1977
#' 
#' @param Z - matrix of integers. Columns are factors, rows are subjects, and integer represent a level for the corresponding factor/column. For a factor with k levels, the corresponding column should take integer values between 1 and k. We assume all levels of all factors appear at least once in the observed data.
#' @param N - number of treatment groups
#' @param DFunc: function to assess variation 
#' @param  weights: vector of factor weights when some factors are more important than others
#' @param formula: choices are a,b,c, corresponding to randomization probability formulas in section 3.3 of paper
#' @param bias: the constant associated with the randomization probability formula in paper. larger values yield less stochasticity. For example, for formula='a', range from 1/N (every assignment equal) to 1 (assignment deterministic). See paper for details. 
#' @param seed: used for assigning treatments and breaking ties
#' @return array X with X[i,j,k] the number of subjects with level j of factor i in treatment k. X[i,j,k]=NA if j>numLevels(factor i).
# assignment vector with assignment[s] the treatment assignment of subject s
# error code: 1 means error, 0 means no errors detected
balanceRandomize = function(Z, N, DFunc = function(x) diff(range(x)), weights=1, formula='a', bias = .7, seed=42) {
  numSubj = nrow(Z) #number of subjects
  M = ncol(Z) #number of factors
  f = apply(Z, 2, max) #number of levels for each factor
  X = array(NA, dim = c(M,max(f),N)) #X[i,j,k] gives number of subjects with level j of factor i in group k
  for (i in 1:M) X[i,1:f[i],] = 0 #NA for j>f[i], 0 otherwise
  dk = vector("double", M) 
  G = vector("double", N) 
  assignment = vector("double", numSubj) 
  p = vector('double', N) #treatment assignment probabilities for subject s
  error = 0
  set.seed(seed)
  for (s in 1:numSubj) { #assign next subject s, given current count array X
    r = Z[s,] #factor levels for subject s
    
    #calculate Gk for each k=1,...,N (see paper)
    for (k in 1:N) {
      
      #calculate dik for each i=1,...,M (see paper)
      for (i in 1:M) {
        Xk = X 
        Xk[i, r[i], k] = Xk[i, r[i], k] + 1 
        dk[i] = DFunc(Xk[i,r[i],])
      }
      
      G[k] = sum(weights*dk)
    }
    
    #    set.seed(42)
    if (formula=='a') {
      bestAssign = which(G==min(G))
      prefer = ifelse(length(bestAssign)==1, bestAssign, sample(bestAssign, 1)) #if ties, break randomly
      p[prefer] = bias
      p[-c(prefer)] = (1-bias)/(N-1) 
    } else if (formula=='c') {
      p = (1/(N-bias))*(1-bias*G/sum(G))
    } else if (formula=='b') {
      p = bias-2*(N*bias-1)/(N*(N+1))*rank(G, ties.method = 'random')
    }
    if (sum(p)!=1) error=1
    assignment[s] = sample(x = 1:N, size = 1, prob = p) #treatment assignment for subject s
    for (i in 1:M) X[i,r[i],assignment[s]] = X[i,r[i],assignment[s]]+1 #update X
  }
  return(list(assignment = assignment, X = X, error = error))
}

  
