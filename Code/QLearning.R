# Implementation of Q-learning to address limitations of the
# qLearn function from the DynTxRegime package namely
# 1) Accommodate combination treatment arms
# 2) Respect infeasible treatment sequences
# Author: John Sperger

#' Q-learning implementation
#' 
#' 
#' 
#' Note: model.type = 'lasso' uses cv to find the tuning parameter value which minimizes the MSE
#' and uses that for predictions
#' 
#' @param indf data frame with study data. Currently requires treatment indicators already be present
#' in the input data (e.g. A2_1 )
#' @param stage1.formula Formula for the stage 1 model. The left hand side of the formula must be PseudoY
#' @param stage2.formula Y ~ X_1*(A2_1 + A2_2),
#' @param stage1.arm.map
#' @param stage2.arm.map
#' @param model.type String which specifies what class of model to use for the mean models. 
#' Possible options are 'lm', 'lasso', and 'rf' 
#' 
#' @param impermissible.arms.df
#' 
#' @return a list where the first two elements are the stage one and stage two models
#' and the third element is a string with the class of model used for fitting the mean models
QLearning <- function(indf,
                      stage1.formula = PseudoY ~ X_1*(A1_1 + A1_2),
                      stage2.formula = Y ~ X_1*(A2_1 + A2_2),
                      stage1.arm.map,
                      stage2.arm.map,
                      model.type = "lm",
                      impermissible.arms.df,
                      ...){
  # Error Checking -------------------------------------------------------------
  
  if (attr(terms((stage1.formula)), which = 'response') == 0) stop("stage1.formula must be a two-sided formula with Pseudo on the left hand side")
  if (attr(terms((stage2.formula)), which = 'response') == 0) stop("stage2.formula must be a two-sided formula with Pseudo on the left hand side")
  if (all.vars(stage1.formula)[1] != "PseudoY") stop("The left hand side of stage1.formula must be 'PseudoY'")
  if (all.vars(stage2.formula)[1] != "Y") stop("The left hand side of stage2.formula must be 'Y'")
  
  model_options <- c("lm", "lasso", "rf")
  
  # Perform some basic pre-processing of the model.type string
  model.type <- tolower(model.type) %>% trimws(.)
  
  if (model.type %in% model_options == FALSE) stop("model.type must be one of 'lm', 'lasso', or 'rf'")
  
  # Modeling -------------------------------------------------------------------
  
  # --- Stage Two
  
  # Fit the second stage model
  
  stage2_data <- .QLearningConstructDataList(in.formula = stage2.formula,
                                             in.data = indf, 
                                             ...)
  
  stage2_fit <- switch(model.type,
                       lm = .QLearningFitLM(data.list = stage2_data),
                       lasso = .QLearningFitLasso(data.list = stage2_data),
                       rf = .QLearningFitRF(data.list = stage2_data))
  
  # Append the counterfactual predictions for stage 2 to the study data
  # Find the pseudovalue to use for stage one regression (taking into account feasibility of arm sequences)
  
  indf_w_preds <- .QLearningCreateCounterfactuals(base.data = indf,
                                          fitted.model = stage2_fit,
                                          model.formula = stage2.formula,
                                          treatment.arm.df = stage2.arm.map) %>% 
    .QLearningFindPseudoVal(indf.with.pred.pos = .,
                            impermissible.trt.pairs = impermissible.arms.df)
    
  # --- Stage One
  
  stage1_data <- .QLearningConstructDataList(in.formula = stage1.formula,
                                             in.data = indf_w_preds, 
                                             ...)
  
  stage1_fit <- switch(model.type,
                       lm = .QLearningFitLM(data.list = stage1_data),
                       lasso = .QLearningFitLasso(data.list = stage1_data),
                       rf = .QLearningFitRF(data.list = stage1_data))
  
  to_return <- list(Stage1Mod = stage1_fit,
                    Stage2Mod = stage2_fit,
                    model.class = class(stage1_fit))
  
  return(to_return)
}


#' This is a function which takes advantage of the long (one row per ptid and arm) 
#' out-of-sample
#' 
#' @param q.mod the resulting object from the \code{QLearning} function which is a list
#' where the models for both stages are present
#' @param oos.data out-of-sample data created by \code{make1OutofSampleSimDataSet}
#'  this will be a long data frame with a row for each out-of-sample patient and hypothetical treatment sequence
#'  
#' @return an {n out-of-sample} by 3 data set where the columns are the ID (ptid),
#' the estimated optimal arm at stage 1 (A1hat) and stage 2 (A2hat)
PredictOptSeqQLearningOOS <- function(q.mod, oos.data) {
  est_seq_df <- oos.data %>% 
    mutate(Q1pred = predict(q.mod$Stage1Mod, oosData),
           Q2pred = predict(q.mod$Stage2Mod, oosData)) %>% 
  group_by(ptid) %>% 
    summarise(A1Hat = A1[which.max(Q1pred)],
              A2Hat = A2[which.max(Q2pred)],
              MuAhat = Mu[which.max(Q1pred + Q2pred)],
              maxMu = maxMu[1],
              DeltaMuAhat = DeltaMu[which.max(Q1pred + Q2pred)],
              PercOracle = max(0, MuAhat/maxMu),
              .groups = "drop_last")
  
  return(est_seq_df)
}

#' This is a function which takes advantage of the long (one row per ptid and arm) 
#' out-of-sample
#' 
#' @param q.mod the resulting object from the \code{QLearning} function which is a list
#' where the models for both stages are present
#' @param oos.data out-of-sample data created by \code{make1OutofSampleSimDataSet}
#'  this will be a long data frame with a row for each out-of-sample patient and hypothetical treatment sequence
#'  
#' @return a single row data frame with columns for Oracle performance measures
PercOracleQLearningOOS <- function(q.mod, oos.data){
  perf_summary <- PredictOptSeqQLearningOOS(q.mod = q.mod, oos.data = oos.data) %>% 
    summarize(MeanVal = mean(MuAhat),
              OracleVal = mean(maxMu),
              PercOracle = mean(MuAhat)/mean(maxMu),
              DifOracle = OracleVal - MeanVal)
  
  return(perf_summary)
}

################################################################################
### Q-learning Core Functions
#
################################################################################

#' Find the pseudovalue to use instead of the observed Y1 using the counterfactual predictions
#' from the stage 2 model and information about the impermissible treatment arm sequences
#' to ensure only valid sequences are considered. 
#' @param indf.with.pred.pos study data with the predicted outcomes under all A2 arm assignments 
#' @param impermissible.trt.pairs data frame where rows are \emph{invalid} treatment sequences and columns are stages (A1, A2)
#' 
#' @return input data frame with two columns appended containing the pseudo value, 'PseudoY', and the estimated
#' optimal stage 2 treatment for that individual 'A2OptHat'
.QLearningFindPseudoVal <- function(indf.with.pred.pos,
                                    impermissible.trt.pairs){
  
  trt_po_names_and_position <- indf.with.pred.pos %>% 
    select(starts_with("YpredA2")) %>% 
    names(.) %>% 
    map_chr(., ~str_replace(., "YpredA2", "")) %>% 
    tibble(A2 = .) %>% mutate(Position = 1:n())
  
  # Reduce the impermissible.trt.pairs df to a single row per A1 arm (or less, if tere are  no impermissible pairs)
  # where 
  impermissible_list_form <- impermissible.trt.pairs %>% 
    left_join(., trt_po_names_and_position, by = "A2") %>% 
    group_by(A1) %>% 
    summarise(ImpermissibleArmNames = list(A2),
              ImpermissibleArmPositions = list(Position),
              .groups = "drop_last")
  
  indf_with_pseudo <- indf.with.pred.pos %>% 
    left_join(., impermissible_list_form, by = "A1") %>% 
    nest(., preds = starts_with("YpredA2")) %>% 
    rowwise %>% 
    mutate(PossiblePreds = list(unlist(preds)[-ImpermissibleArmPositions]),
           A2OptPos = which.max(unlist(PossiblePreds)),
           A2OptHat = str_replace(names(PossiblePreds)[A2OptPos], "YpredA2", ""),
           PseudoY = max(unlist(PossiblePreds))) %>% 
    unnest(cols = preds) %>% 
    select(-PossiblePreds, -A2OptPos)
  
  return(indf_with_pseudo)
  }

#' Create a matrix of counterfactual predictions for a given stage. 
#' 
#' @param fitted.model
#' @param base.data
#' @param model.formula
#' @param treatment.arm.df
#' 
.QLearningCreateCounterfactuals <- function(base.data,
                                            fitted.model,
                                            model.formula,
                                            treatment.arm.df){
  
  stopifnot(CheckArmsMatchDataArmDF(study.data = base.data, treatment.arm.map = treatment.arm.df) == TRUE)
  
  # Get the arm variable name which should be either A1 or A2. Will be used for naming columns of the counterfactuals
  arm_var_name <- treatment.arm.df %>% select(-contains("_")) %>% names
  stopifnot(length(arm_var_name) == 1 & arm_var_name %in% c("A1", "A2"))
  
  fitted_model_type <- class(fitted.model)
  
  data_arms_removed <- base.data %>% select(-any_of(colnames(treatment.arm.df)))
  
  # Pre-allocate placeholder matrix for speed
  prediction_mat <- matrix(nrow = nrow(base.data), ncol = nrow(treatment.arm.df))
  colnames(prediction_mat)  <- paste0("Temp", 1:ncol(prediction_mat))
  
  # Create predictions under each treatment assignment
  for(cur_index in 1:nrow(treatment.arm.df)){
    counterfactual_data <- bind_cols(data_arms_removed, treatment.arm.df[cur_index,]) %>% 
      .QLearningConstructCovariateMatrix(in.formula = model.formula,
                                  in.data = .)
    
    prediction_mat[, cur_index] <- switch(fitted_model_type,
      lm = predict(fitted.model, newdata = as.data.frame(counterfactual_data)),
      cv.glmnet = glmnet:::predict.cv.glmnet(fitted.model, 
                                         newx = as.matrix(counterfactual_data[, -str_detect(colnames(counterfactual_data), "(Intercept)")]), #For whatever reason glmnet doesn't want x to have an intercept column
                                         s = "lambda.min"),
      randomForest = randomForest:::predict.randomForest(fitted.model, newdata = counterfactual_data))
    
    cur_arm_name <- treatment.arm.df %>% slice(cur_index) %>%  pull(all_of(arm_var_name))
    
    colnames(prediction_mat)[cur_index] <- paste0("Ypred", arm_var_name, cur_arm_name)
  }
  
  # The column names of the predictions get messed up unless you first turn the matrix into a data frame or tibble
  data_w_preds <- prediction_mat %>% as_tibble(.) %>% bind_cols(base.data, .)
  
  return(data_w_preds)
}


################################################################################
### Q-learning Mean Modeling Functions
#
################################################################################

#'  Fit a (non-regularized) linear model to the data
#'  @param data.list a named two-element list as created by \code{.QLearningConstructDataList}
#'  where the first element, called 'x', is the covariate matrix
#'  and the second element 'y' is the corresponding response vector. 
#'  
#'  @return an lm object
.QLearningFitLM <- function(data.list, ...){
  
  # This way of calling lm preserves the covariate names
  lm_fit <- lm(data.list$y ~ . + 0, 
               data = as.data.frame(data.list$x),
               ...)
  
  return(lm_fit)
}

#'  Fit an L1 penalized linear model to the data
#'  @param data.list a named two-element list as created by \code{.QLearningConstructDataList}
#'  where the first element, called 'x', is the covariate matrix
#'  and the second element 'y' is the corresponding response vector. 
#'  
#'  @return a fitted glmnet object
.QLearningFitLasso <- function(data.list, ...){
  
  glmnet_fit <- glmnet::cv.glmnet(x = data.list$x,
                           y = data.list$y,
                           ...)
  
  return(glmnet_fit)
}

#'  Fit a random forest to the data
#'  @param data.list a named two-element list as created by \code{.QLearningConstructDataList}
#'  where the first element, called 'x', is the covariate matrix
#'  and the second element 'y' is the corresponding response vector. 
#'  
#'  @return a fitted randomForest object
.QLearningFitRF <- function(data.list, ...){
  
  rf_fit <- randomForest::randomForest(x = data.list$x,
                               y = data.list$y,
                               ...)
  
  return(rf_fit)
}

################################################################################
### Q-learning Helper Functions - Covariate Construction
#
# Mean modeling functions use a model matrix and response vector input format
################################################################################

#' This function takes \code{glmnetUtil}'s \code{makeModelComponents} function
#' and modifies it to remove redundant factor levels if desired.
#' Returns a list with the data matrix and the response vector
#' @param in.formula
#' @param in.data
#' @param names.to.remove
#' @param create.intercept
#' 
#' @return a two-element list where the first element is the covariate matrix and the second element is the response vector
.QLearningConstructDataList <- function(in.formula, 
                                        in.data, 
                                        names.to.remove = NULL,
                                        create.intercept = FALSE,
                                   ...){
  
  model.comp <- glmnetUtils:::makeModelComponents(formula = in.formula,
                                                  data = in.data)
  
  x_mat <- model.comp$x[, !colnames(model.comp$x) %in% names.to.remove]
  
  # If there is only a single row model.comp$x is a vector not a matrix
  if (nrow(in.data) == 1) {
    x_mat <- matrix(data = x_mat, nrow = 1)
    colnames(x_mat) <- colnames(model.comp$x)[!colnames(model.comp$x) %in% names.to.remove]
  } 
  
  y_vec <- model.comp$y
  
  if (create.intercept == TRUE){
    x_mat <- cbind(rep(1, times = nrow(x_mat)), x_mat)
    colnames(x_mat)[1] <- "(Intercept)"
  } 
  
  if (any(duplicated(colnames(x_mat))) == TRUE) stop("Formula resulted in duplicate column names, check that there aren't repeated terms in the formula")
  
  lasso_data_list <- list(x = x_mat,
                          y = y_vec)
  
  return(lasso_data_list)
}

#' This function takes \code{glmnetUtil}'s \code{makeModelComponents} function
#' and modifies it to allow users to specify variables to remove. Examples of variables 
#' include redundant factor levels.
.QLearningConstructCovariateMatrix <- function(in.formula, 
                                               in.data, 
                                               names.to.remove = NULL,
                                               create.intercept = TRUE, ...){
  return(.QLearningConstructDataList(in.formula = in.formula,
                                in.data = in.data,
                                names.to.remove = names.to.remove,
                                create.intercept = create.intercept)$x)
}
