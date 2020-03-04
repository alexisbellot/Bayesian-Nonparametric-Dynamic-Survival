# A Bayesian Approach to Modelling Longitudinal Data in Electronic Health Records

This is an R implementation of the paper ["A Bayesian Approach to Modelling Longitudinal Data in Electronic Health Records"](https://arxiv.org/abs/1912.09086). 

In this project we present a new approach to model time-to-event data (such as death or infection times in a medical context) using repeated measurements of covariates recorded over time. Our objective is to integrate the heterogeneity in medical data, by which we mean including irregularly sampled measurements, informatively missing data and demographic information among others, to improve predictions of survival. In our paper we discuss the relevance of these issues. In this repository we develop a Bayesian nonparametric model that aims to model all this information, quantify the uncertainty around model predictions in a principled manner at the individual level, while avoiding to make assumptions on the data generating process and adapting model complexity to the structure of the data. We believe both contributions to be particularly important for personalizing health care decisions as predictions may be uncertain due to lack of data, while we can expect the underlying heterogeneous physiological time series to vary wildly across patients.

*Please cite the above paper if this resource is used in any publication.*

## Requirements

* R version 3.5 or later.
* Packages: "pec","bartmachine", "dynpred", "JM".

## First steps
To get started, check Demo.R which will guide you through an application of our algorithm with publicly available data, which we have preprocessed and can be downloaded from our data folder. 

If you have questions or comments about anything regarding this work, please do not hesitate to contact [Alexis](https://alexisbellot.github.io/Website/)
