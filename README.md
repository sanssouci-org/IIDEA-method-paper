# IIDEA-method-paper : Powerful and interpretable control of false discoveries in differential expression studies

This repository contains scripts to reproduce all experiments and figures of the IIDEA method paper (https://hal.archives-ouvertes.fr/hal-03601095/). 

# How to reproduces figures 

Each figures of the paper are in the folder "figures". To reproduce each figures, use the script as follow : 

| Figures  | Scripts         | Options |
| :--------------- |:---------------| :-----|
| Figure 01   |   ce texte  |           |
| Figure 02  | figure-02.R  |        |
| Figure 03  | figure-03_S-01.R          |    |
| Figure 04  | 00_utils.R > 01_setup.R > 02_run_differntial_expression.R > 03_plot_JER-control-results.R >  | technology = "RNAseq" |
| Figure S-01  | figure-03_S-01.R |    Aligné à droite |
| Figure S-02  | comparison_tests.R |      |
| Figure S-03  | permutation_studies.R          |    technology = "RNAseq" |
| Figure S-04  | permutation_studies.R          |    technology = "RNAseq" |
| Figure S-05  | 00_utils.R > 01_setup.R > 02_run_JER-control-results.R | technology = "RNAseq" |
| Figure S-06  | sample_studies_RNAseq.R        |         |
| Figure S-06  | sample_studies_RNAseq.R        |         |
| Figure S-05  | 00_utils.R > 01_setup.R > 02_run_JER-control-results.R | technology = "microarray" |
| Figure S-06  | sample_studies_microarray.R        |         |
| Figure S-06  | continuous_covariate.R        |         |
