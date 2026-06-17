README - 10_Machinelearning/code

10_check_RF_precision_metrics_MeJAvsControl_SAvsControl.Rmd: loads Boruta and Cross-validation results to extract RF-important ASVs, make plots on RF ASVs individual importance (justification of Ds-top quartile cut-off), make plots on RF precision and kappa.

Machine_Learning_custom_functions.R: several custom functions used on the machine learning pipeline. most importantly:
    * they facilitate running Boruta as a feature selection algorithm to classify samples according to a metadata column (Such as "Stress")
    * they put the whole spliting, training and testing pipeline with caret on a single function. this process is replicated to accound for random smapling when spiting training and testing datasets
