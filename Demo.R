###
### Example performance estimation on the publicly available PBC2 data 
###

## Load utility functions

source('/utils.R')


## Minimal example of survival estimates and plot

data = get_pbc_data(landmark_time=2,horizon=5)
n = nrow(data)
predictions = survival_estimates(LMdata_train = data[1:(n/2),],LMdata_test=data[(n/2):n,])
survival = predictions[[1]] # survival predictions
times = predictions[[3]] # corresponding times
plot(times,survival[1,],type='l') # survival probabilities over time for patient 1

## Run algorithm and compute performance

perf = cross_validated_performance(LMpoints=c(2),horizon=5,k=2)



