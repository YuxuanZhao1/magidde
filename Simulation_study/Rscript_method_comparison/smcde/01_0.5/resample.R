check.weights = function(weights,log=FALSE,normalized=TRUE) 
{
  stopifnot(length(weights)>0)
  
  if (all(log,normalized))
    warning("logged and normalized weights are unusual")
  
  if (!normalized) 
  {
    weights = renormalize(weights, log, engine="R") # returns unlogged
    log = FALSE
  }
  
  if (log) 
  {
    weights=exp(weights)
  } else
  {
    if (any(weights<0)) 
    {
      warning("log=FALSE, but negative there is at least one negative weight.\nAssuming log=TRUE.")
    }
  }
  
  return(weights)  
}

inverse.cdf.weights = function(weights, uniforms=runif(length(weights)), engine="R")
{
  check.weights(weights)
  check.weights(uniforms)
  
  engine=pmatch(engine, c("R","C"))
  if (is.unsorted(uniforms)) uniforms = sort(uniforms)
  n.samples = length(uniforms)
  
  switch(engine,
         {
           # R implementation
           ids       = integer(n.samples)
           cusum     = cumsum(weights)
           index     = 1
           for (i in 1:n.samples) 
           {
             found = FALSE
             while (!found) 
             {
               if (uniforms[i] > cusum[index]) 
               {
                 index = index + 1
               }
               else 
               {
                 found = TRUE
               }
             }
             ids[i] = index
           }
           return(ids)
         },
         {
           # C implementation
           out = .C("inverse_cdf_weights_R", 
                    as.integer(length(weights)),
                    as.double(weights),
                    as.integer(n.samples),
                    as.double(uniforms), 
                    id=integer(n.samples))
           return(out$id)
         })  
}

systematic.resample = function(weights, num.samples=length(weights), engine="R")
{
  check.weights(weights, log=F, normalized=T)
  stopifnot(num.samples>0)
  
  engine=pmatch(engine, c("R","C"))
  n = length(weights)
  
  switch(engine,
         {
           u = runif(1,0,1/num.samples)+seq(0,by=1/num.samples,length=num.samples)
           return(inverse.cdf.weights(weights,u,engine="R"))
         },
         {
           # C implementation
           out = .C("systematic_resample_R", 
                    as.integer(n),
                    as.double(weights),
                    as.integer(num.samples),
                    id = integer(num.samples))
           return(out$id+1)
         })
}


