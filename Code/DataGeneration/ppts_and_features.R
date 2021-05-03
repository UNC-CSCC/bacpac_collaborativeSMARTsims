# Functions for making participants and their features 
# Author: Nikki Freeman
# Last modified on 22 November 2020

# Make participants -----------------------------------------------------------
#' Make study participants
#'
#' @param N integer; total number of study participants
#'
#' @return dataframe with patient identifiers
makePpts <- function(N){
  out <- data.frame(ptid = 1:N)
  
  return(out)
}

# Make features/covariates ----------------------------------------------------


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
covariateFn_Mvn <- function(df, numBinaryCovars, props,
                            mu.vec, sigma.mat, d.partial.obs, prop.unobserved){
  out <- df %>% 
    makeBinaryCovars(numBinaryCovars, props) %>%
    MakeMVNCovariates(mu.vec, sigma.mat, d.partial.obs, prop.unobserved)
  
  return(out)
}

#' Add binary covariates to the generated patients that take values in {-1, 1}
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
        mutate(!!sym(covarNames[i]) := 2*rbinom(n = 1, size = 1, prob = props[i])-1) %>%
        ungroup()
    }
    return(df)
  }
}

#' Add binary covariates to the generated patients that takes value in {0, 1}
#'
#' @param df dataframe; at a minimum contains the same number of rows 
#' as participants
#' @param numCovars integer; number of binary covariates to add to df
#' @param props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#'
#' @return dataframe with the generated binary covariates appended
makeBinaryCovarsZeroOne <- function(df, numCovars, props){
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


# Make features/covariates ----------------------------------------------------
# Partially-observed covariates -----------------------------------------------
# -----------------------------------------------------------------------------

#' Add multivariate normal covariates to the generated ppts
#' Optionally have some of the covaraites be partially observed. 
#' Creates indicators of whether the covariate was observed for the partially observed covariates
#' @param df dataframe; at a minimum contains the same number of rows as
#' participants 
#' @param mu.vec vector; means of normal covariates
#' @param sd sigma.mat; sds of normal covariates
#'
#' @return dataframe; df with the normal covariates appended. 
#' W_x for continuous covariates
#' U_x for partially observed covariates
#' O_x for indicator of whether the subject had expensive covariates measured
MakeMVNCovariates <- function(df, mu.vec, sigma.mat, d.partial.obs = 0, prop.unobserved){
  if (is.null(mu.vec)) return(df)
  
  cont_covar <- mvtnorm::rmvnorm(n = nrow(df), mean = mu.vec, sigma = sigma.mat)
  
  if (d.partial.obs > 0){
    d.observed <- ncol(cont_covar) - d.partial.obs
    d.true.params <- ncol(cont_covar)
    
    # Want to keep the partially observed variables for data generation
    obs_var_names <- paste0("W_", 1:length(mu.vec))
    partial_var_names <- paste0("U_", 1:d.partial.obs)
    obs_indicator_var_names <- paste0("O_", 1:d.partial.obs)
    
    start_pos_partial_covar <- min(d.true.params, d.true.params - d.partial.obs + 1)
    
    partial_covar <- matrix(cont_covar[, start_pos_partial_covar:d.true.params], ncol = d.partial.obs)
    partial_covar[sample(1:nrow(partial_covar), 
                         size = floor((prop.unobserved)*nrow(partial_covar))), ] <- NA

    obs_indicators <- is.na(partial_covar) == FALSE
      
    cont_covar <- cbind(cont_covar, partial_covar, obs_indicators)
    
    
    colnames(cont_covar) <- c(obs_var_names, partial_var_names, obs_indicator_var_names)
    
    
  } else{ colnames(cont_covar) <- paste0("W_", 1:ncol(cont_covar))}
  
  df <- bind_cols(df, as_tibble(cont_covar))
  
  return(df)
}


# Make features/covariates ----------------------------------------------------
# Cluster-correlated data generation-------------------------------------------
# -----------------------------------------------------------------------------



#' Make covariate function that make
#'
#' @param df dataframe; each row should correpsond to a ppt 
#' @param n.sites number of sites
#' @param numBinaryCovars integer; non-neg number of binary covariates to add
#' @param props vector of probabilities of successes for each binary covariate
#' @param numNormalCovars integer; non-neg number of normal covariates to add
#' @param mu vector of means for the normal covariates
#' @param sd vector of sds for the normal covariates
#'
#' @return dataframe; df with the covariates appended
#' 
#' @export
covariateFn_Cluster <- function(df, n.sites, 
                                bin.props, bin.iccs,
                                mu.vec, overall.sd, normal.iccs, ...){
  out <- df %>% 
    MakeClusters(n.sites = n.sites) %>% 
    MakeCorrelatedBinaryCovars(bin.props, bin.iccs, ...) %>%
    MakeClusterCorrelatedNormalCovars(mu.vec, overall.sd, normal.iccs)
  
  return(out)
}


#' Create clusters (enrollment sites)
#' @param df the current data frame
#' @param n.sites the number of sites that will be enrolling patients (the number of clusters)
#' 
#' @return the input data frame with a column indicating site membership appended
MakeClusters <- function(df, n.sites){
  return(df %>% mutate(Site = sample(1:n.sites, size = n(), replace = TRUE)))
}

#' Create cluster-correlated binary covariates using the \code{fabricatr} package. 
#' By default covariates use a {-1, 1} encoding
#' @param df the current data frame
#' @param bin.props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#' @param bin.iccs vector of intra-cluster correlation coefficients the same length as bin.props, or a single constant that will be used for all covariates
#' @param zero.one.encoding logical indicating whether {0,1} encoding should be used instead of the default {-1, 1} encoding. Defaults is FALSE
#' 
#' @return the input data frame with the covariates appended
MakeCorrelatedBinaryCovars <- function(df, bin.props, bin.iccs, zero.one.coding = FALSE){
  stopifnot("Site" %in% colnames(df))
  
  num_covars <- length(bin.props)
  
  if(num_covars == 0) return(df)
  
  # Recycle bin.iccs so its length matches num_covars
  if(length(bin.iccs) < num_covars){
    if(length(bin.iccs) > 1) warning("Length of bin.iccs and bin.props do not match, recycling bin.iccs. Bin.iccs is not a constant, verify this mismatch is intended")
    bin.iccs <- rep(bin.iccs, length.out = num_covars)
  }
  
  # Create covariates
  bin_covar_mat <- matrix(nrow = nrow(df), ncol = num_covars)
  
  for(i in 1:num_covars){
    bin_covar_mat[,i] <- fabricatr::draw_binary_icc(prob = bin.props[i],
                                                    clusters = df$Site,
                                                    ICC = bin.iccs[i])
  }
  
  # Transform to {-1, 1} coding unless zero.one.coding is specified
  if(zero.one.coding == FALSE) bin_covar_mat <- 2*bin_covar_mat - 1
  
  # Naming
  covar_names <- paste("X", 1:num_covars, sep = "_")
  colnames(bin_covar_mat) <- covar_names
  
  df_with_bin_covars <- bind_cols(df, as_tibble(bin_covar_mat))
  
  return(df_with_bin_covars)
}


#' Create cluster-correlated binary covariates using the \code{fabricatr} package. 
#' By default covariates use a {-1, 1} encoding
#' @param df the current data frame
#' @param mu.vec vector of means. The length of this vector determines the number of normal covariates to generate
#' @param overall.sd the standard deviation for each covariate. The within and between cluster variance is derived from this and the ICC 
#' @param normal.iccs vector of intra-class correlation coefficients
#' 
#' @return the input data frame with the covariates appended
MakeClusterCorrelatedNormalCovars <- function(df, mu.vec, overall.sd, normal.iccs){
  stopifnot("Site" %in% colnames(df))
  
  num_covars <- length(mu.vec)
  
  if(num_covars == 0) return(df)
  
  # Recycle normal.iccs so its length matches num_covars
  if(length(normal.iccs) < num_covars){
    if(length(normal.iccs) > 1) warning("Length of bin.iccs and bin.props do not match, recycling bin.iccs. Bin.iccs is not a constant, verify this mismatch is intended")
    normal.iccs <- rep(normal.iccs, length.out = num_covars)
  }
  
  # Create covariates
  norm_covar_mat <- matrix(nrow = nrow(df), ncol = num_covars)
  
  for(i in 1:num_covars){
    norm_covar_mat[,i] <- fabricatr::draw_normal_icc(mean = mu.vec[i],
                                                     clusters = df$Site,
                                                     ICC = normal.iccs[i],
                                                     total_sd = overall.sd)
  }
  
  # Naming
  covar_names <- paste("W", 1:num_covars, sep = "_")
  colnames(norm_covar_mat) <- covar_names
  
  df_with_norm_covars <- bind_cols(df, as_tibble(norm_covar_mat))
  
  return(df_with_norm_covars)
}
