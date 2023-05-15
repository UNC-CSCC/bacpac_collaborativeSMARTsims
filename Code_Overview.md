BACPAC Code Overview
========================================================
author: 
date: 2021-02-15
autosize: true

Overview
========================================================

Code Sections

- Data Generation
- Out-of-sample Data
- Analysis

Each section is controlled by two lists, `metadata` and `args` which dictate how everything works.

Someday this may move to a JSON file, for now ¯\\_(ツ)_/¯

The code is not at all optimized. Still be mindful when writing code, but the aim was flexibility.

Data Generation Overview
========================================================


```r
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
```

Data Generation Overview - Metadata
========================================================


```r
metadata <- list(N = 1000,
                 firstLineTreatments = as.character(0:3),
                 secondLineTreatments = as.character(0:3),
                 augmentationTreatments = as.character(4:9),
                 standardOfCareTreatment = "0")

impermissible_arm_seqs <- bind_rows(expand_grid(A1 = c("1", "2", "3"), A2 = "0"),
                                    expand_grid(A1 = "0", A2 = c(as.character(1:6))),
                                    expand_grid(A1 = "1", A2 = c("6", "8", "9")),
                                    expand_grid(A1 = "2", A2 = c("5", "7", "9")),
                                    expand_grid(A1 = "3", A2 = c("4", "7", "8")))
          
metadata$possibleTreatmentSequences <- generatePossibleTreatmentSequencesGrid(
  first_line_trts = metadata$firstLineTreatments,
  second_line_trts = metadata$secondLineTreatments,
  augment_trts = metadata$augmentationTreatments,
  impermissible_trt_pairs_df = impermissible_arm_seqs)
```


Data Generation Overview - Metadata
========================================================


```r
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
```

Common "reserved" names
========================================================
| Name      | Used for |
| ----------- | ----------- |
| A1      | character denoting stage-one treatment assignment   |
| A2   | character denoting stage-two treatment assignment     |
| A1_#   | 0/1 indicator if treatment # was received at stage 1. Not necessarily the same as A1 because there are combination treatments |
| Y1     | Stage-one outcome   |
| Y2   | Stage-two outcome    |
| Y      | Overall outcome (relationship to Y1 and Y2 depends on simulation setting)  |


Combination Treatments
========================================================

Stage Two example


```
 A2 A2_0 A2_1 A2_2 A2_3
  0    1    0    0    0
  1    0    1    0    0
  2    0    0    1    0
  3    0    0    0    1
  4    0    1    1    0
  5    0    1    0    1
  6    0    0    1    1
  7    1    1    0    0
  8    1    0    1    0
  9    1    0    0    1
```


Contributing to the codebase
========================================================
 Use [Google's R Style Guide](https://google.github.io/styleguide/Rguide.html)

- snake_case for variable names 
- BigCamelCase for functions

Qualify the namespace for all functions except for `tidyverse` functions and base R functions


```r
# This
foo <- glmnet::glmnet(bar)
  
# Not this
foo <- glmnet(bar)
```
