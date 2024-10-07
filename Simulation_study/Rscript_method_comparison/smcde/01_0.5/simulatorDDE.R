## Simulator ###
#set.seed(452)
#library(deSolve)

parameters <- c(r = 0.8, P = 2, tau = 3)
parameters

initial <- c(3500)
initial
#state
DE_model <- function(t, y, parameters) {
  ## assign parameter values to parameter variables
  ## the function 'unname()' just removes the name of the parameter - this is unnecessary, it just cleans things up a bit
  r <- unname(parameters['r'])
  P <- unname(parameters['P'])
  tau <- unname(parameters['tau'])
  
  ## define time lags ###
  if (t < tau)
    lag <- initial
  else
    lag <- lagvalue(t - tau)
  
  dy <- r*y*(1-lag/(1000*P))
  
  return(list(dy, dy = dy))

}


## Simulate the system!
output <- dede(y=initial, times=times, func=DE_model, parms=parameters, atol = 1e-10)


## Simulate observations y ##
y <- exp(as.vector(log(output[,2]) + rnorm(length(times), 0, sigma1)))

