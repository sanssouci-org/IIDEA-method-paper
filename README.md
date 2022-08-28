# Powerful and interpretable control of false discoveries in differential expression studies

This repository contains scripts to reproduce all the numerical experiments and figures of the paper:

> Enjalbert-Courrech, N., & Neuvial, P. (2022). *Powerful and interpretable control of false discoveries in differential expression studies*. https://hal.archives-ouvertes.fr/hal-03601095/ 

A pdf version of each figure is stored in "figures/", and the R code to reproduce the figures is found in "scripts/".

Note that some of the numerical experiments take a long time to run since we are performing 1000 simulations based on real data sets with 10,000 to 20,000 genes by 100 to 200 observations. The experiments are parallelized using the `future` R package. The number of 'workers' used for parallel computing can be set as follows before running a given experiment:

```r
library(future)
plan(multisession, workers = 2) # parallelize on 2 nodes
```
## Calibration algorithm

- Figure 1

  - `source("scripts/algo-calibration.R")`

## RNAseq data analysis

- `source("scripts/perf_RNAseq_run.R")` to produce the following figures:

  - Figures 2 and S-1 (confidence curves)
  - Figure 3 (volcano plot)
  - Figure S-2 (limma vs Wilcoxon p-values)

## Performance evaluation: JER control and power

- Figures 4 and S-5 (RNAseq data)

  - `source("scripts/perf_RNAseq_run.R")` to run the experiments
  - `source("scripts/perf_RNAseq_plot.R")` to produce the figures

- Figures S-8 (microarray data)

  - `source("scripts/perf_microaray_run.R")` to run the experiments
  - `source("scripts/perf_microarray.R")` to produce the figure

## Influence of sample size

- Figures S-6 and S-7 (RNAseq data)

  - `source("scripts/sample-size_RNAseq_run.R")` to run the experiments
  - `source("scripts/sample-size_RNAseq_plot.R")` to produce the figures

- Figures S-9 (microarray data)

  - `source("scripts/sample-size_microaray_run.R")` to run the experiments
  - `source("scripts/sample-size_microarray.R")` to produce the figure

## Influence of the number of permutations

- Figures S-3 and S-4

  - `source("scripts/number-of-permutations_run.R")` to run the experiments
  - `source("scripts/number-of-permutations_plot.R")` to produce the figures

## Continuous covariates

- Figures S-10

  - `source("scripts/continuous-covariate_run.R")` to run the experiments
  - `source("scripts/continuous-covariate_plot.R")` to produce the figure
