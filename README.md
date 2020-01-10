
---
title: Improving disaggregation models of malaria incidence by ensembling non-linear models of prevalence.
---


This is the code used to run the models and create figures for the manuscript `Improving disaggregation models of malaria incidence by ensembling non-linear models of prevalence`.
Unfortunately, due to data sharing agreements we cannot distribute the data so this code is not reproducible.


The aim of this study is to test one method for improving disaggregation regression models with a focus on malaria.
Here, we use malaria prevalence surveys to run an initial suite of machine learning models.
Predictions from those models are used as covariates in the disaggregation regression model.
This two stage model is compared to a baseline model in which environmental covariates are used directly.


The scripts `run_machine_learn.R` and `run_machine_learn_*data.R` run the machine learning models and need to be run first.

The the scripts `run_analysis_*.R` run all the different disaggregation models, and cross-validation runs, for each country.

Finally, `make_figures_and_summaries.R` creates the figures.

`joint_model.cpp` is the file that defines the disaggregation model and is used by TMB to fit the models.



