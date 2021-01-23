library(tidyverse)
library(DynTxRegime)

library(tictoc)
library(future)
library(furrr)
library(here)

#if (str_detect(here::here(), "jsperger") == TRUE) plan(multicore, workers = 8)
plan(sequential)
scripts <- str_subset(list.files(path = "../Code", full.names = TRUE), ".R$")
quietly(map(scripts, source))


################################################################################
## Settings
##
################################################################################


##------------------------------------------------------------------------------
## Settings: Simulation Run meta-settings
##------------------------------------------------------------------------------
# Seed is from the NC Pick 4 Daytime Draw Sunday, Jan 17  https://nclottery.com/Pick4
seed_to_use <-2814
set.seed(seed_to_use)


# Whether to generate new simulated study data or read in existing 
# If False, it will try to read in an .Rds of existing simulated data from read_sim_file_string
run_sim_flag <- TRUE
read_sim_file_string <- NULL

if (run_sim_flag == TRUE & is_null(read_sim_file_string) == FALSE) stop("run_sim_flag must be FALSE if read_sim_file_string is specified")
if (run_sim_flag == FALSE & is_null(read_sim_file_string) == TRUE) stop("read_sim_file_string must be specified if run_sim_flag is FALSE")

# Whether the simulated study data should be saved (to an .Rds)
save_sim_data_flag <- FALSE
save_sim_file_string <- paste0("../SimRuns/best_guess_no_subgroups_seed",seed_to_use, ".RDS")

if (save_sim_data_flag == TRUE & is_null(save_sim_file_string)) stop("Must specify save_sim_file_string if save_sim_data_flag is TRUE")

# Whether to generate out-of-sample data 
# If read_sim_file_string is specified, it will try to read in an .Rds of existing out-of-sample
# simulated data from read_oos_file_string
gen_oos_flag <- TRUE

read_oos_file_string <- NULL

if (gen_oos_flag == TRUE & is_null(read_oos_file_string) == FALSE) stop("gen_oos_flag must be FALSE if read_oos_file_string is specified")


# Whether the simulated out-of-sample data should be saved (to an .Rds)
save_oos_data_flag <- FALSE
save_oos_file_string <- "../SimRuns/OOS_best_guess_no_subgroups.RDS"

if (save_sim_data_flag == TRUE & is_null(save_sim_file_string)) stop("Must specify save_sim_file_string if save_sim_data_flag is TRUE")


# Whether to fit models to and analyze the simulated study data
# If , it will try to read in an .Rds of existing simulated data from the same scenario
run_analysis_flag <- TRUE
read_analysis_file_string <- NULL

if (run_analysis_flag == TRUE & is_null(read_analysis_file_string) == FALSE) stop("run_analysis_flag must be FALSE if read_analysis_file_string is specified")

save_analysis_flag <- FALSE
save_analysis_file_string <-  "../SimRuns/performance_best_guess_no_subgroups.RDS"

# Create L data sets per setting
L <- 1000

##------------------------------------------------------------------------------
## Settings: Simulated Study Data Generation
##------------------------------------------------------------------------------
metadata <- list(N = 1000,
                 firstLineTreatments = as.character(0:3),
                 secondLineTreatments = as.character(0:3),
                 augmentationTreatments = as.character(4:6),
                 standardOfCareTreatment = "0")

impermissible_arm_seqs <- bind_rows(expand_grid(A1 = c("1", "2", "3"), A2 = "0"),
                                    expand_grid(A1 = "0", A2 = c("4", "5", "6")),
                                    tibble(A1 = c("1", "2", "3"), 
                                           A2 = c("6", "5", "4")))
          
metadata$possibleTreatmentSequences <- generatePossibleTreatmentSequencesGrid(
  first_line_trts = metadata$firstLineTreatments,
  second_line_trts = metadata$secondLineTreatments,
  augment_trts = metadata$augmentationTreatments,
  impermissible_trt_pairs_df = impermissible_arm_seqs)

# Read in parameter vectors from CSVs
stage1_param_df <- read.csv("./ParameterCSVs/SubgroupSc1-ATEscenario1/sg1_stage1.csv")

stage1_param_vec <- stage1_param_df$Coefficient
names(stage1_param_vec) <- stage1_param_df$Parameter

stage2_param_df <- read.csv("./ParameterCSVs/SubgroupSc1-ATEscenario1/sg1_stage2.csv")

stage2_param_vec <- stage2_param_df$Coefficient
names(stage2_param_vec) <- stage2_param_df$Parameter


#W_1 is being used as a random intercept

stage1_args <- list( stage.suffix = 1,
                     treatment.arm.df = stage_1_trts <- tibble(A1 = as.character(0:3),
                                                               A1_1 = c(0, 1, 0, 0),
                                                               A1_2 = c(0, 0, 1, 0),
                                                               A1_3 = c(0, 0, 0, 1)),
                     arm.var.name = "A1",
                     true.model.formula = formula(~ 1 + W_1 + A1_1 + A1_2 + A1_3 + A1_2:X_1), 
                     true.param.vec = stage1_param_vec,
                     sigma.noise = sqrt(.5),
                     include.trt.indicators = TRUE,
                     include.mu = FALSE)

stage2_args <- list(stage.suffix = 2,
                    treatment.arm.df = tibble(A2 = as.character(0:6),
                                              A2_1 = c(0, 1, 0, 0, 1, 1, 0),
                                              A2_2 = c(0, 0, 1, 0, 1, 0, 1),
                                              A2_3 = c(0, 0, 0, 1, 0, 1, 1)),
                    arm.var.name = "A2",
                    true.model.formula = formula(~ 1 + W_1 + A2_1 + A2_2 + A2_3 + A2_1:A2_2 + A2_1:A2_3 + A2_2:A2_3 + A2_2:X_1), 
                    true.param.vec = stage2_param_vec,
                    sigma.noise = sqrt(.5),
                    include.trt.indicators = TRUE,
                    include.mu = FALSE)

args <- list(
  makePpts_fn = "makePpts",
  makePpts_args = list(metadata[["N"]]),
  makeCovariates_fn = "covariateFn_v1",
  makeCovariates_args = list(numBinaryCovars = 3, props = c(0.5, 0.2, 0.4),
                             numNormalCovars = 1, mu = 0, sd = sqrt(.5)),
  allocateStage1Treatments_fn = "StratifiedBlockRandomizationStage1",
  allocateStage1Treatments_args = list(first.line.trts = metadata$firstLineTreatments,
                                       strata.vars.syms = c(sym("X_1"))),
  generateY1_fn = "GenNormalOutcomeByFormula",
  generateY1_args = stage1_args,
  assignResponderStatus_fn = "assignResponderStatusByQuantile",
  assignResponderStatus_args = list(ctsOutcome1 = "Y1", 
                                    cutoffs = c(0.5, 0.75, .875),
                                    status.labels = c("bad", "medium", "good", "excellent")),
  allocateStage2Treatments_fn = "allocationFn_stage2_trt_grid_4resps",
  allocateStage2Treatments_args = list(trt.options.grid = metadata$possibleTreatmentSequences,
                                       lowest_augmentation_trt_number = min(as.numeric(metadata$augmentationTreatments))),
  generateY2_fn = "GenNormalOutcomeByFormula",
  generateY2_args = stage2_args,
  generateY_fn = "addY1AndY2",
  generateY_args = list(NA))


##------------------------------------------------------------------------------
## Settings: Scenario Modifications
##------------------------------------------------------------------------------

metadata_settings_scenario1 <- map(c(800, 1000, 1200), ~list_modify(metadata, N = .))

metadata_list <- c(metadata_settings_scenario1)

args_list_scenario1 <- map(c(800, 1000, 1200), ~list_modify(args, makePpts_args = list(.)))

args_list <- c(args_list_scenario1)

##------------------------------------------------------------------------------
## Settings: Out-of-sample Data Generation
##------------------------------------------------------------------------------

oos_metadata <-  list(nOOS = 1e3,
                      firstLineTreatments = metadata$firstLineTreatments,
                      secondLineTreatments = metadata$secondLineTreatments,
                      augmentationTreatments = metadata$augmentationTreatments,
                      standardOfCareTreatment = metadata$standardOfCareTreatment)

oos_args <- list(
  makePpts_fn = args$makePpts_fn,
  makePpts_args = list(N = oos_metadata[["nOOS"]]),
  makeCovariates_fn = args$makeCovariates_fn,
  makeCovariates_args = args$makeCovariates_args,
  MakePOs_fn = "calcTrueMeansAndPOs",
  makePOs_args = list(trt_options_grid = metadata$possibleTreatmentSequences,
                      outcome_functions_list = list(args$generateY1_fn, 
                                                    args$generateY2_fn,
                                                    args$assignResponderStatus_fn,
                                                    args$generateY_fn),
                      outcome_args_list = list(args$generateY1_args,
                                               args$generateY2_args,
                                               args$assignResponderStatus_args,
                                               args$generateY_args)))

##------------------------------------------------------------------------------
## Settings: Data Analysis
##------------------------------------------------------------------------------

analysis_args <- list(AnalysisModeling_fn = "QLearning",
                      AnalysisModeling_args = list(stage1.formula = formula(PseudoY ~ 1 + A1_1 + A1_2 + A1_3 + A1_2:X_1),
                                                   stage2.formula = formula(Y ~ 1 + I(A1_1 + A2_1) + I(A1_2 + A2_2) + I(A1_3 + A2_3) + 
                                                                              I(A1_2 + A2_2):X_1 + 
                                                                              A2_1:A2_2 + A2_1:A2_3 + A2_2:A2_3),
                                                   stage1.arm.map = args$generateY1_args$treatment.arm.df,
                                                   stage2.arm.map = args$generateY2_args$treatment.arm.df,
                                                   model.type = "lm",
                                                   impermissible.arms.df = impermissible_arm_seqs),
                      PredictSequence_fn = "predictTreatSequence",
                      PredictSequence_args = NULL,
                      makeInSampleMetrics_fns = NULL,
                      makeInSampleMetrics_args_list = NULL,
                      makeOOSMetrics_fns = "",
                      makeOOSMetrics_args = list())

################################################################################
## Data Generation
##
################################################################################

##------------------------------------------------------------------------------
## Single study example for debugging
##------------------------------------------------------------------------------

if (FALSE) {
  s1 <- makeSimDataWrapper(metadata = metadata, args = args, 1) %>% 
    mutate(Trt1 = A1_1 + A2_1, 
           Trt2 = A1_2 + A2_2, 
           Trt3 = A1_3 + A2_3)
  
  m1 <- lm(Y ~ 1 +  A1_1 + A1_2 + A1_3 + A2_1 + A2_2 + A2_3 + A2_1:A2_2 + A2_1:A2_3 + A2_2:A2_3 + A1_2:X_1 + A2_2:X_1, data = s1)
  m1.simp <- lm(Y ~ Trt1 + Trt2 + Trt3 + A2_1:A2_2 + A2_1:A2_3 + A2_2:A2_3 + Trt2:X_1, data = s1)
  m1.2 <- lm(Y2 ~ 1 + A2_1 + A2_2 + A2_3 + A2_1:A2_2 + A2_1:A2_3 + A2_2:A2_3 + A2_2:X_1, data = s1)
  m1.1 <- lm(Y1 ~ 1 + A1_1 + A1_2 + A1_3 + A1_2:X_1 , data = s1)
  
  q1 <- exec(analysis_args$AnalysisModeling_fn, indf = s1, !!!analysis_args$AnalysisModeling_args)
  
}
##------------------------------------------------------------------------------
## Data Generation: Out-of-sample data
##------------------------------------------------------------------------------

if (gen_oos_flag == TRUE) {
  print("Generating out-of-sample data")
  
  tic()
  # This is going to be a long data frame because there are {#unique treatment sequences} rows per participant
  oosData <-  make1OutofSampleSimDataSet(metadata = oos_metadata,
                                         args = oos_args) %>% 
    left_join(., args$generateY1_args$treatment.arm.df, by = "A1") %>% 
    left_join(., args$generateY2_args$treatment.arm.df, by = "A2") %>% 
    mutate(`A1_2:X_1` = A1_2*X_1,
           `A2_2:X_1` = A2_2*X_1)
  toc()
}

if (is_null(read_oos_file_string) == FALSE) oosData <- readRDS(read_oos_file_string)

if (save_oos_data_flag == TRUE) saveRDS(oosData, file = save_oos_file_string)

##------------------------------------------------------------------------------
## Data Generation: Run the simulation
##------------------------------------------------------------------------------

if (run_sim_flag == TRUE){
  print("Generating simulated datasets")
  tic()
  sim_data_list <-   furrr::future_map2(.x = metadata_list, 
                                        .y = args_list,
                                        ~map(1:L, makeSimDataWrapper, metadata = .x, args = .y),
                                        .options = furrr::furrr_options(seed = TRUE))
  toc()
} else{ sim_data_list <- readRDS(read_sim_file_string)}

if (save_sim_data_flag == TRUE) saveRDS(sim_data_list, save_sim_file_string)

if (run_analysis_flag == TRUE) {
  if (FALSE) {
    
  
  print("Q-learning")
  tic()
  q_mods_list <-   map(sim_data_list, 
                       ~furrr::future_map(.x = .,
                                     ~exec(analysis_args$AnalysisModeling_fn, indf = ., !!!analysis_args$AnalysisModeling_args),
                                     .options = furrr::furrr_options(seed = TRUE)))
  toc()
  
  tic()
  
  ### Warning - hard coding here
  # TODO - Remove hard coding
  oracle_summary <- map(q_mods_list, 
                        ~map_dfr(.x = ., ~PercOracleQLearningOOS(q.mod = ., oos.data = oosData))) %>% 
    bind_rows(., .id = "Scenario") %>% 
    left_join(., y = tibble(Scenario = as.character(1:3), N = c(800, 1000, 1200)), by = "Scenario")
  
  
  toc()
  }

  
  print("Hypothesis testing")


  tic()
  stage1_power_results_sc1 <- map(sim_data_list[1:length(metadata_settings_scenario1)],
                                  ~map_dfr(.x = ., ~CalcStage1WaldUnadjustedPvals(study.data = .,
                                                                                  resp.formula = formula(Y1 ~ A1),
                                                                                  coefs.to.test = c("A11", "A12", "A13")))) %>% 
    map(., ~AdjustPvals(., adj.method = "BH"))
  
  stage1_pval_summary <- map_dfr(stage1_power_results_sc1, ~colMeans(. < .05))
  
  toc()
  
 #analysis_results <- list(oracle_summary, stage1_power_results_sc1)
  
  # 
} 

if (is_null(read_analysis_file_string) == FALSE) analysis_results <- readRDS(read_analysis_file_string)

if (save_analysis_flag == TRUE) saveRDS(analysis_results, save_analysis_file_string)

plan(sequential)


##------------------------------------------------------------------------------
## Scratch work
##------------------------------------------------------------------------------


if (FALSE){ 
  
  
  oos_preds <- map(q_mods, ~PredictOptSeqQLearningOOS(q.mod = ., oos.data = oosData))
  
  stage1_wald <-map(sim_data_list[[3]], ~CalcStage1WaldUnadjustedPvals(study.data = .,
                                                         resp.formula = formula(Y1 ~ A1_1 + A1_2 + A1_3),
                                                         contrast.mat = diag(4)[-1,]))
  
  CalcStage1WaldUnadjustedPvals <- function(study.data,
                                            resp.formula = formula(Y1 ~ A1*X_1),
                                            contrast.mat = NULL,
                                            coefs.to.test = NULL)
  
  oosData <- oosData %>% 
    mutate(Q1pred = predict(q1$Stage1Mod, newdata = oosData),
           Q2pred = predict(q1$Stage2Mod, newdata = oosData)) %>% 
    select(-.groups)
  
  useful <- oosData %>% select(A1, A2, ptid, Mu, DeltaMu, Q1pred, Q2pred)
  
  dtr_assigned_trts <-useful %>% group_by(ptid) %>% summarise(A1Hat = A1[which.max(Q1pred)],
                                                              A2Hat = A2[which.max(Q2pred)],
                                                              MuAhat = Mu[which.max(Q1pred + Q2pred)],
                                                              DeltaMuAhat = DeltaMu[which.max(Q1pred + Q2pred)],)
  
  rf_args <- list_modify(analysis_args, AnalysisModeling_args = list_modify(analysis_args$AnalysisModeling_args, model.type = "rf"))
  
  
  rf_q_mods_list <-   map(sim_data_list, 
                          ~furrr::future_map(.x = .,
                                             ~exec(rf_args$AnalysisModeling_fn, indf = ., !!!rf_args$AnalysisModeling_args),
                                             .options = furrr::furrr_options(seed = TRUE)))
  
  rf_oracle_summary <- map(rf_q_mods_list, 
                        ~map_dfr(.x = ., ~PercOracleQLearningOOS(q.mod = ., oos.data = oosData))) %>% 
    bind_rows(., .id = "Scenario") %>% 
    left_join(., y = tibble(Scenario = as.character(1:3), N = c(800, 1000, 1200)), by = "Scenario")
}