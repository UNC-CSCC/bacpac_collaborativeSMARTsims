library(tidyverse)
library(DynTxRegime)

library(future)
library(furrr)
library(here)

if (str_detect(here::here(), "jsperger") == TRUE) plan(multicore, workers = 8)

scripts <- str_subset(list.files(path = "../Code", full.names = TRUE), ".R$")
quietly(map(scripts, source))


##------------------------------------------------------------------------------
## Simulation Run Settings
##------------------------------------------------------------------------------
set.seed(8936)
# NC Pick 4 Daytime Draw Monday, Dec 7
# https://nclottery.com/Pick4

# Whether to generate new simulated study data. 
# If False, it will try to read in an .Rds of existing simulated data from the same scenario
run_sim_flag <- TRUE

# Whether the simulated study data should be saved (to and .Rds)
save_sim_data_flag <- TRUE

# Whether to fit models to and analyze the simulated study data
# If False, it will try to read in an .Rds of existing simulated data from the same scenario
run_analysis_flag <- TRUE
save_analysis_flag <- TRUE

# Create L data sets per setting
L <- 2000

##------------------------------------------------------------------------------
## Simulation - Scenario Settings
##------------------------------------------------------------------------------
metadata <- list(N = 600,
                 firstLineTreatments = as.character(0:3),
                 secondLineTreatments = as.character(0:4),
                 standardOfCareTreatment = "0")

metadata$possibleTreatmentSequences <- generatePossibleTreatmentSequencesGrid(
  first_line_trts = metadata$firstLineTreatments,
  second_line_trts = metadata$secondLineTreatments,
  augment_trts = metadata$augmentationTreatments,
  impermissible_trt_pairs_df = expand_grid(A1 = c("1", "2", "3"), A2 = "0"))

args <- list(
  makePpts_fn = "makePpts",
  makePpts_args = list(metadata[["N"]]),
  makeCovariates_fn = "covariateFn_v1",
  makeCovariates_args = list(numBinaryCovars = 3, props = c(0.6, 0.2, 0.4),
                             numNormalCovars = 1, mu = 0, sd = 1),
  allocateStage1Treatments_fn = "allocationFn_stage1_v1",
  allocateStage1Treatments_args = list(metadata$firstLineTreatments),
  generateY1_fn = "generateNormalY1",
  generateY1_args = list(intercept = 0, 
                         coefs = c("A1_1*X_1" = 0.15, 
                                   "A1_2*X_1" = -.15), 
                         sigma = sqrt(0.5)),
  assignResponderStatus_fn = "assignResponderStatusByQuantile",
  assignResponderStatus_args = list(ctsOutcome1 = "Y1", cutoffs = c(0.2, 0.8)),
  allocateStage2Treatments_fn = "allocationFn_stage2_trt_grid",
  allocateStage2Treatments_args = list(trt.options.grid = metadata$possibleTreatmentSequences),
  generateY2_fn = "generateY2Fn_v1",
  generateY2_args = list(coefs = c("A2_3" = 0.15),
                         sigma = sqrt(0.5)),
  generateY_fn = "setYtoY2",
  generateY_args = list(NA))

trts_including_non_feasible <- generatePossibleTreatmentSequencesGrid(
  first_line_trts = metadata$firstLineTreatments,
  second_line_trts = metadata$secondLineTreatments,
  augment_trts = metadata$augmentationTreatments)

## Out-of-sample

oos_metadata <-  list(nOOS = 1e3,
                      firstLineTreatments = metadata$firstLineTreatments,
                      secondLineTreatments = metadata$secondLineTreatments,
                      standardOfCareTreatment = metadata$standardOfCareTreatment)

oos_args <- list(
  makePpts_fn = args$makePpts_fn,
  makePpts_args = list(N = oos_metadata[["nOOS"]]),
  makeCovariates_fn = args$makeCovariates_fn,
  makeCovariates_args = args$makeCovariates_args,
  MakePOs_fn = "calcTrueMeansAndPOs",
  makePOs_args = list(trt_options_grid = trts_including_non_feasible,
                      outcome_functions_list = list(args$generateY1_fn, 
                                                    args$generateY2_fn,
                                                    args$assignResponderStatus_fn,
                                                    args$generateY_fn),
                      outcome_args_list = list(args$generateY1_args,
                                               args$generateY2_args,
                                               args$assignResponderStatus_args,
                                               args$generateY_args)))


# This is going to be a long data frame
oosData <-  make1OutofSampleSimDataSet(metadata = oos_metadata,
                                       args = oos_args)


## Model specification
moMain_stage2 <- modelObj::buildModelObj(model = ~ Y1,
                                         solver.method = 'lm')

moCont_stage2 <- modelObj::buildModelObj(model = ~ 1,
                                         solver.method = 'lm')

moMain_stage1 <- modelObj::buildModelObj(model = ~ -1,
                                         solver.method = 'lm')

moCont_stage1 <- modelObj::buildModelObj(model = ~ X_1,
                                         solver.method = 'lm')

moCont_stage1_scenario2 <- modelObj::buildModelObj(model = ~ 1,
                                                   solver.method = 'lm')
## Analysis
analysis_metadata <-  list(N = metadata$N,
                           nOOS = oos_metadata$nOOS,
                           firstLineTreatments = metadata$firstLineTreatments,
                           secondLineTreatments = metadata$secondLineTreatments,
                           standardOfCareTreatment = metadata$standardOfCareTreatment)

analysis_args <- list(AnalysisModeling_fn = "doQLearningListReturn",
                      AnalysisModeling_args = list("moMain_stage2" = moMain_stage2, 
                                                   "moCont_stage2" = moCont_stage2,
                                                   "moMain_stage1" = moMain_stage1, 
                                                   "moCont_stage1" = moCont_stage1,
                                                   "fSetFunction" = fSetFunScenario3),
                      PredictSequence_fn = "predictTreatSequence",
                      PredictSequence_args = NULL,
                      makeInSampleMetrics_fns = NULL,
                      makeInSampleMetrics_args_list = NULL,
                      makeOOSMetrics_fns = "",
                      makeOOSMetrics_args = list())

## Scenario 2
analysis_args_scenario2 <- list_modify(analysis_args, AnalysisModeling_args = list("moMain_stage2" = moMain_stage2, 
                                                                                   "moCont_stage2" = moCont_stage2,
                                                                                   "moMain_stage1" = moMain_stage1, 
                                                                                   "moCont_stage1" = moCont_stage1_scenario2,
                                                                                   "fSetFunction" = fSetFunScenario3))

generateY1_args_scenario2 = list(intercept = 0, 
                                 coefs = c("A1_1" = 0.15), 
                                 sigma = sqrt(0.5))

generateY2_args_scenario2 = list(coefs = c("A2_1" = 0.15),
                                 sigma = sqrt(0.5))

args_scenario2 <- list_modify(args, generateY1_args = generateY1_args_scenario2,
                              generateY2_args = generateY2_args_scenario2) 

oos_args_scenario2 <- list_modify(oos_args, 
                                  makePOs_args = list(trt_options_grid = trts_including_non_feasible,
                                                      outcome_functions_list = list(args$generateY1_fn, 
                                                                                    args$generateY2_fn,
                                                                                    args$assignResponderStatus_fn,
                                                                                    args$generateY_fn),
                                                      outcome_args_list = list(args_scenario2$generateY1_args,
                                                                               args_scenario2$generateY2_args,
                                                                               args$assignResponderStatus_args,
                                                                               args$generateY_args)))

oosDataScenario2 <-  make1OutofSampleSimDataSet(metadata = oos_metadata,
                                                args = oos_args_scenario2)

##------------------------------------------------------------------------------
## Create Modifications to Vary
##------------------------------------------------------------------------------

metadata_settings_scenario1 <- map(c(600, 800, 1000, 1200), ~list_modify(metadata, N = .))
metadata_settings_scenario2 <- map(c(600, 800, 1000, 1200), ~list_modify(metadata, N = .))

metadata_list <- c(metadata_settings_scenario1,
                   metadata_settings_scenario2)

args_list_scenario1 <- map(c(600, 800, 1000, 1200), ~list_modify(args, makePpts_args = list(.)))
args_list_scenario2 <- map(c(600, 800, 1000, 1200), ~list_modify(args_scenario2, makePpts_args = list(.)))

args_list <- c(args_list_scenario1,
               args_list_scenario2)


##------------------------------------------------------------------------------
## Run the simulation
##------------------------------------------------------------------------------

if (run_sim_flag == TRUE){
  library(tictoc)
  tic()
  sim_data_list <-   furrr::future_map2(.x = metadata_list, 
                                        .y = args_list,
                                        ~map(1:L, makeSimDataWrapper, metadata = .x, args = .y),
                                        .options = furrr::furrr_options(seed = TRUE))
  toc()
} else{ sim_data_list <- readRDS("sim_data_list_design3.RDS")}

if (save_sim_data_flag == TRUE) saveRDS(sim_data_list, "sim_data_list_design3.RDS")

if (run_analysis_flag == TRUE) {
  
  tic()
  analysis_results_sc1 <- map(sim_data_list[1:length(metadata_settings_scenario1)],
                              ~analyzeSimulationRunsPercOracleWrapper(metadata = analysis_metadata, 
                                                                      args = analysis_args,
                                                                      study.data.list = .,
                                                                      oos.data = oosData)) %>% 
    bind_rows(.) %>% mutate(Design = "3", Scenario = "1")
  toc()
  
  tic()
  analysis_results_sc2 <- map(sim_data_list[(length(metadata_settings_scenario1)+1):length(metadata_list)],
                              ~analyzeSimulationRunsPercOracleWrapper(metadata = analysis_metadata, 
                                                                      args = analysis_args,
                                                                      study.data.list = .,
                                                                      oos.data = oosDataScenario2)) %>% 
    bind_rows(.) %>% mutate(Design = "3", Scenario = "2")
  
  toc()
  
  analysis_results <- bind_rows(analysis_results_sc1, analysis_results_sc2)
  
  
} else{
  analysis_results <- readRDS("design3_analysis.rds")
}
if (save_analysis_flag == TRUE) saveRDS(analysis_results, "design3_analysis.rds")

plan(sequential)
