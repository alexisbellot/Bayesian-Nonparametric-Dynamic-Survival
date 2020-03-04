# This script provides utility functions to be used for dynamic survival prediction

# Required packages
options(java.parameters="-Xmx10000m")
library(dynpred)
library(bartMachine)
library(pec)
library(JM)
library(ggplot2)
set_bart_machine_num_cores(2)


## Function to construct longitudinal data -----------------------------------------------------------------------

make_LM_data <- function(data, outcome, LM, horizon, covs,
                         format = c("wide","long"), id, rtime, right=TRUE){
  if (missing(id))
    stop("argument 'id' should be specified for long format data")
  if (missing(rtime))
    stop("argument 'rtime' should be specified for long format data")
  ord <- order(data[[id]],data[[rtime]])
  data <- data[ord,]
  ids <- unique(data[[id]])
  n <- length(ids)
  # initialize LMdata; copy first row of each subject
  LMdata <- data[which(!duplicated(data[[id]])),]
  for (i in 1:n) {
    wh <- which(data[[id]]==ids[i])
    di <- data[wh,]
    #print(i)
    idx <- cut(LM,c(data[[rtime]][wh],Inf),right=right,labels=FALSE)
    if (!is.na(idx)) LMdata[i,] <- di[idx,]
    else {
      LMdata[i,] <- di[1,]
      LMdata[[covs$varying]][i] <- NA
      LMdata[[rtime]][i] <- NA
    } 
  }
  
  LMdata <- LMdata[LMdata[[outcome$time]] > LM,] # keep patients that suervive beyond LM
  if (format=="long") LMdata <- LMdata[!is.na(LMdata[[id]]),]
  # apply administrative censoring at horizon
  LMdata[outcome$status] <- LMdata[[outcome$status]] *
    as.numeric(LMdata[[outcome$time]] <= horizon)
  LMdata[outcome$time] <- pmin(as.vector(LMdata[[outcome$time]]),horizon)
  LMdata$LM <- LM
  if (format=="long")
    cols <- match(c(id,outcome$time,outcome$status,covs$fixed,covs$varying,rtime,"LM"),
                  names(LMdata))
  else
    cols <- match(c(outcome$time,outcome$status,covs$fixed,covs$varying,"LM"),
                  names(LMdata))
  return(LMdata[,cols])
}

## Load pbc data -----------------------------------------------------------------------------

get_pbc_data = function(landmark_time, horizon){
  
  # import data (pre-processed)
  pbc_data = read.csv(file='C:/Users/abellot/Documents/R Scripts/Dynamic_Survival_Experiments/Toy example/pbc.csv')
  
  # convert all variables to numeric
  
  LMdata <- make_LM_data(data=pbc_data,outcome=list(time="Survival",status="Status"),
                         LM=landmark_time,horizon=landmark_time+horizon,
                         covs=list(fixed=c("id","drug","age","sex"),
                                   varying=c("ascites","hepatomegaly", "spiders","noEdema","noDiu",
                                             "serBilir","serChol","albumin","alkaline","SGOT","platelets",
                                             "prothrombin","histologic")),
                         format="long", id="id", rtime="follow_up",right=FALSE)
  LMdata$Survival     = round(LMdata$Survival, digits = 1)
  LMdata              = subset(LMdata, select=-c(id,id.1,follow_up,LM))
  return(LMdata)
}

## Function to prepare data for BART --------------------------------------------------------------

surv.pre.bart <- function(
  times,
  ## vector of survival times
  
  delta,
  ## vector of event indicators: 1 event, 0 censoring
  
  x.train=NULL,
  ## matrix of covariate regressors
  ## can be NULL, i.e. KM analog
  
  x.test=NULL
  ## matrix of covariate regressors at X.test settings
  
) {
  ## currently does not handle time dependent Xs
  ## can be extended later
  ## most likely via the alternative counting process notation
  
  N <- length(times)
  
  events <- unique(sort(times))
  ## time grid of events including censoring times
  
  K <- length(events)
  
  y.train <- integer(N) ## y.train is at least N long
  
  k <- 1
  
  for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
    y.train[k] <- delta[i]*(times[i] == events[j])
    
    k <- k+1
  }
  
  m <- length(y.train)
  
  
  
  
  if(length(x.train)==0) {
    p <- 0
    n <- 1
    
    X.train <- matrix(nrow=m, ncol=1, dimnames=list(NULL, 't'))
  } else {
    p <- ncol(x.train)
    
    if(length(x.test)>0) n <- nrow(x.test)
    
    X.train <- matrix(nrow=m, ncol=p+1)
    
    if(length(dimnames(x.train)[[2]])>0)
      dimnames(X.train)[[2]] <- c('t', dimnames(x.train)[[2]])
    else dimnames(X.train)[[2]] <- c('t', paste0('x', 1:p))
  }
  
  k <- 1
  
  for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
    if(p==0) X.train[k, ] <- c(events[j])
    else X.train[k, ] <- c(events[j], x.train[i, ])
    
    k <- k+1
  }
  
  if(p==0 | length(x.test)>0) {
    X.test <- matrix(nrow=K*n, ncol=p+1, dimnames=dimnames(X.train))
    
    for(i in 1:n) for(j in 1:K) {
      if(p==0) X.test[j, ] <- c(events[j])
      else X.test[(i-1)*K+j, ] <- c(events[j], x.test[i, ])
    }
  }
  else X.test <- NULL
  
  return(list(y.train=y.train, X.train=X.train, X.test=X.test, times=events, K=K))
}

## Mean survival estimates BART from longitudinal data -----------------------------------------------

survival_estimates = function(LMdata_train,LMdata_test,confidence=F,partial_dependence=F,
                              pd_variable="ascites",prior_proposals=NULL){
  # Compute posterior survival curves with BNDS
  # Inputs:
  # - LMdata_train: output of "make_LM_data", longitudinal data up to a certain time point
  # - LMdata_test: testing data similar to LMdata_train
  pre                         = surv.pre.bart(times=LMdata_train$Survival,delta=LMdata_train$Status,
                                              x.train= as.matrix(subset(LMdata_train,select=-c(Survival,Status))),
                                              x.test=as.matrix(subset(LMdata_test,select=-c(Survival,Status))))
  X.train                     = data.frame(pre$X.train)
  X.test                      = data.frame(pre$X.test)
  
  bart                        = bartMachine(X=X.train, y=as.factor(pre$y.train), num_iterations_after_burn_in = 1000,
                                            use_missing_data = TRUE,verbose=F,num_trees = 50, cov_prior_vec = prior_proposals)
  
  inclusion = get_var_props_over_chain(bart)
  posterior_proposals     = as.numeric(get_var_props_over_chain(bart))
  
  cat("BART done ", "\n")
  
  if(confidence){
    # credible intervals, numeric matrix with 2 rows and nrow(X.test) columns
    posterior                   = apply(bart_machine_get_posterior(bart, X.test)$y_hat_posterior_samples,1,
                                        function(x)quantile(x,probs=c(0.05,0.95)))}
  
  
  if(partial_dependence == TRUE){
    pd = pd_plot(bart, pd_variable, prop_data = 0.5)
    pd = pnorm(pd$bart_avg_predictions_by_quantile[10])-pnorm(pd$bart_avg_predictions_by_quantile[1])}
  else pd = NULL
  
  # Arraging output for survival
  bart$surv.test              = predict(bart, X.test)
  if(confidence){
  bart$lower                  = posterior[1,] #lower bound of cred. interval
  bart$upper                  = posterior[2,]} #upper bound of cred. interval
  unique.times                = unique(sort(pre$X.train[ , 1]))
  K                           = length(unique.times)
  H                           = nrow(pre$X.test)/K
  for(h in 1:H){ for(j in 2:K){
    bart$surv.test[K*(h-1)+j] = bart$surv.test[K*(h-1)+j-1]*bart$surv.test[K*(h-1)+j]
    
    if(confidence){
      bart$upper[K*(h-1)+j]     = bart$upper[K*(h-1)+j-1]*bart$upper[K*(h-1)+j]
      bart$lower[K*(h-1)+j]     = bart$lower[K*(h-1)+j-1]*bart$lower[K*(h-1)+j]}}}
  
  if(confidence){
    upper                       = t(matrix(bart$upper, ncol=H))
    lower                       = t(matrix(bart$lower, ncol=H))}
  else upper = lower = NULL
  
  survbartmachine             = t(matrix(bart$surv.test, ncol=H))
  
  return(list(survbartmachine,inclusion,unique.times,upper,lower,pd))
}

## Mean survival estimates from fitted model -----------------------------------------------

create_model_data = function(LMdata_train,LMdata_test){
  # Create data for BNDS
  # Input:
  # - LMdata: output of make_LM_data
  pre                         = surv.pre.bart(times=LMdata_train$Survival,delta=LMdata_train$Status,
                                              x.train= as.matrix(subset(LMdata_train,select=-c(Survival,Status))),
                                              x.test=as.matrix(subset(LMdata_test,select=-c(Survival,Status))))
  X.train                     = data.frame(pre$X.train)
  X.test                      = data.frame(pre$X.test)
  unique.times                = unique(sort(pre$X.train[ , 1]))
  
  return(X.test, unique.times)
}

## Mean survival estimates from fitted model -----------------------------------------------

predict_survival_estimates = function(model, LMdata,confidence=F,partial_dependence=F,
                              pd_variable="ascites",prior_proposals=NULL){
  # Compute posterior survival curves with BNDS
  # Inputs:
  # - LMdata: output of "create_model_data", longitudinal data in a format amenable to BNDS
  X.test = LMdata[1]
  unique.times = LMdata[2]
  
  inclusion = get_var_props_over_chain(model)
  posterior_proposals     = as.numeric(get_var_props_over_chain(model))
  
  cat("BART done ", "\n")
  
  if(confidence){
    # credible intervals, numeric matrix with 2 rows and nrow(X.test) columns
    posterior                   = apply(bart_machine_get_posterior(model, X.test)$y_hat_posterior_samples,1,
                                        function(x)quantile(x,probs=c(0.05,0.95)))}
  
  
  if(partial_dependence == TRUE){
    pd = pd_plot(model, pd_variable, prop_data = 0.5)
    pd = pnorm(pd$bart_avg_predictions_by_quantile[10])-pnorm(pd$bart_avg_predictions_by_quantile[1])}
  else pd = NULL
  
  # Arraging output for survival
  model$surv.test              = predict(bart, model)
  if(confidence){
    bart$lower                  = posterior[1,] #lower bound of cred. interval
    bart$upper                  = posterior[2,]} #upper bound of cred. interval
  #unique.times                = unique(sort(pre$X.train[ , 1]))
  K                           = length(unique.times)
  H                           = nrow(X.test)/K
  for(h in 1:H){ for(j in 2:K){
    model$surv.test[K*(h-1)+j] = model$surv.test[K*(h-1)+j-1]*model$surv.test[K*(h-1)+j]
    
    if(confidence){
      model$upper[K*(h-1)+j]     = model$upper[K*(h-1)+j-1]*model$upper[K*(h-1)+j]
      model$lower[K*(h-1)+j]     = model$lower[K*(h-1)+j-1]*model$lower[K*(h-1)+j]}}}
  
  if(confidence){
    upper                       = t(matrix(model$upper, ncol=H))
    lower                       = t(matrix(model$lower, ncol=H))}
  else upper = lower = NULL
  
  survbartmachine             = t(matrix(model$surv.test, ncol=H))
  
  return(list(survbartmachine,inclusion,unique.times,upper,lower,pd))
}


## Performance estimates for BART and Cox -----------------------------------------------------------------------

cross_validated_performance = function(LMpoints, horizon, k=2){
  
  # require surv.pre.bart function
  # data must be one-hot-oncoding format, all variables numeric
  # LM points are the time points at which predictions are desired
  # horizon is maximum prediction horizon
  
  cindex_cox = brier_cox = sd.cindex_cox = sd.brier_cox = list()
  cindex_bart = brier_bart = sd.cindex_bart = sd.brier_bart = list()
  priors = NULL
  
  for( i in 1:length(LMpoints)){
    
    # Create LM data set
    LMdata = get_pbc_data(LMpoints[i],horizon)
    n                   = nrow(LMdata)
    folds               = split(sample(n), seq_len(k))
    
    
    xval.fold = function(fold,prior_update=F) {
      
      # BART performance
      predictions = survival_estimates(LMdata_train = LMdata[-fold,],LMdata_test=LMdata[fold,],
                                       prior_proposals = priors)
      
      survbartmachine     = predictions[[1]]
      #priors              = predictions[[2]] + 1e-8
      unique.times        = predictions[[3]]
      if(prior_update)
        return(priors)
      
      else
        cox = fit_cox(LMdata_train = LMdata[-fold,])
      nonzero.coef = cox[[2]]
      cox          = cox[[1]]
      
      ## Performance metrics
      cindex_fold_cox = pec::cindex(object=cox,formula=Surv(Survival,Status)~1,data=LMdata[fold,c(TRUE,TRUE,nonzero.coef)], 
                                    eval.times=unique.times)$AppCindex
      brier_fold_cox  = pec::pec(object=cox,formula=Surv(Survival,Status)~1,data=LMdata[fold,], start=unique.times[1],
                                 maxtime=unique.times[length(unique.times)],times = unique.times,exact=FALSE,reference=FALSE)$AppErr
      cindex_fold_bart = pec::cindex(object=survbartmachine,formula=Surv(Survival,Status)~1,data=LMdata[fold,], 
                                     eval.times=unique.times)$AppCindex
      brier_fold_bart  = pec::pec(object=survbartmachine,formula=Surv(Survival,Status)~1,data=LMdata[fold,], start=unique.times[1],
                                  maxtime=unique.times[length(unique.times)],times = unique.times,exact=FALSE,reference=FALSE)$AppErr
      
      return(data.frame('cindex_fold_cox'=cindex_fold_cox,'brier_fold_cox'=brier_fold_cox,
                        'cindex_fold_bart'=cindex_fold_bart,'brier_fold_bart'=brier_fold_bart,
                        'unique.times'=unique.times))
    }
    
    #priors = xval.fold(folds[[1]],prior_update = T)
    results               = lapply(folds, xval.fold)
    results               = do.call("rbind", results)
    
    
    cindex_cox[[i]]           = tapply(results[,1],as.factor(results[,5]),mean)
    brier_cox[[i]]            = tapply(results[,2],as.factor(results[,5]),mean)
    sd.cindex_cox[[i]]        = tapply(results[,1],as.factor(results[,5]),sd)
    sd.brier_cox[[i]]         = tapply(results[,2],as.factor(results[,5]),sd)
    cindex_bart[[i]]           = tapply(results[,3],as.factor(results[,5]),mean)
    brier_bart[[i]]            = tapply(results[,4],as.factor(results[,5]),mean)
    sd.cindex_bart[[i]]        = tapply(results[,3],as.factor(results[,5]),sd)
    sd.brier_bart[[i]]         = tapply(results[,4],as.factor(results[,5]),sd)
    names(cindex_cox[[i]])    = names(brier_cox[[i]]) = names(sd.cindex_cox[[i]]) = names(sd.brier_cox[[i]]) = unique(sort(results[,5]))
    names(cindex_bart[[i]])    = names(brier_bart[[i]]) = names(sd.cindex_bart[[i]]) = names(sd.brier_bart[[i]]) = unique(sort(results[,5]))
    
  }
  return(list('cindex_cox'=cindex_cox, 'brier_cox'=brier_cox, 'sd.cindex_cox'=sd.cindex_cox,'sd.brier_cox'=sd.brier_cox,
              'cindex_bart'=cindex_bart, 'brier_bart'=brier_bart, 'sd.cindex_bart'=sd.cindex_bart,'sd.brier_bart'=sd.brier_bart))
}

## Fit Cox ---------------------------------------------------------------------------------------------

fit_cox = function(LMdata_train){
  
  glmnet.cv <- cv.glmnet(as.matrix(subset(LMdata_train,select=-c(Survival,Status))),
                         Surv(LMdata_train$Survival,LMdata_train$Status),family="cox")
  glmnet.obj <- glmnet.cv$glmnet.fit
  optimal.lambda <- glmnet.cv$lambda.min    # For a more parsimoneous
  lambda.index <- which(glmnet.obj$lambda==optimal.lambda) 
  # take beta for optimal lambda 
  optimal.beta  <- glmnet.obj$beta[,lambda.index] 
  # find non zero beta coef 
  nonzero.coef <- abs(optimal.beta)>0 
  selectedVar   <- LMdata_train[,c(TRUE,TRUE,nonzero.coef)] 
  cox = coxph(Surv(Survival,Status) ~.,data=selectedVar,x=TRUE)
  
  if(any(is.na(cox$coefficients)))  cox$coefficients[is.na(cox$coefficients)]=0
  
  return(list(cox,nonzero.coef))
}



## Get variable importance over time ------------------------------------------------------------------

variable_inclusion = function(LMpoints,horizon){
  
  inclusion = list()
  
  for(i in 1:length(LMpoints)){
    
    LMdata = get_CF_data(LMpoints[i],horizon)
    predictions = survival_estimates(LMdata_train = LMdata,LMdata_test=LMdata[1:10,],
                                     prior_proposals = NULL)
    inclusion[[i]]      = predictions[[2]]
  }
  inclusion = data.frame(t(data.frame(inclusion)))
  rownames(inclusion) = LMpoints
  return(inclusion)
}

