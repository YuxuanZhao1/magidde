### adaptive annealing SMC, use B-spline to solve ODE###
ESSVec <- list()

parametersStore <- list()
sigmaStore <- list()
#initialStore <- list()
tempSeq <- list()
tempDifSeq <- list()
cStore <- list()
lambdaStore <- list()


logZ = 0
alpha <- 0
alphaDif <- 0

particles <- cbind(c, r_particles, P_particles, tau_particles, sigma, lambda)

W <- rep(1/NP, NP)
logW <- log(W)
iter <- 0

while(alpha < 1){
  iter = iter + 1
  print(iter)
  c <- particles[,1:ncol(basismat)] 
  sigma <- particles[,ncol(basismat)+4] 
  parameters <- particles[,(ncol(basismat)+1):(ncol(basismat)+3)]
  logRatio <- rep(NA, NP)
  for(i in 1:NP){
    logRatio[i] <- tempPosterior(as.vector(parameters[i,]), sigma[i], y, as.vector(c[i,]), lambda[i])
  }
  alphaDif <- bisection(func, 0, 1, W, logRatio, CESSthresholds)
  
  alpha <- alpha + alphaDif
  print(alpha)
  if(alpha > 1){
    alpha <- 1
    alphaDif <- 1- tempSeq[[iter-1]]
  }
  tempSeq[[iter]] <- alpha
  tempDifSeq[[iter]] <- alphaDif
  
  logw <- alphaDif*logRatio + logW
  logmax = max(logw)
  W = exp(logw-logmax)/sum(exp(logw-logmax)) 
  logW = log(W)
  logZ = logZ + log(sum(exp(logw - logmax))) + logmax 
  print(logZ)
  RESS <- ESS(logW)
  ESSVec[[iter]] <- RESS

  if(RESS < resampleThreshold){
    ancestor <- systematic.resample(W, num.samples = length(W), engine="R")
    particesOld = particles[ancestor,]
    logW = rep(-log(NP), NP)
    W = rep(1/NP, NP)
  }else{
    particesOld = particles
  }
  
  ptm <- proc.time()
  
  c <- particles[, 1:ncol(basismat)]
  parameters <- particles[, (ncol(basismat)+1):(ncol(basismat)+3)]
  sigma <- particles[, (ncol(basismat)+4)]
  
  for(i in 1:NP){
    ## Evaluate Basis functions ###
    Basis_eval <- newBasis(knots, parameters[i,3], bsbasis)
    quadpts_eval <- Basis_eval$quadpts
    quadwts_eval <- Basis_eval$quadwts
    D0quadbasismat_eval <- Basis_eval$D0quadbasismat
    D1quadbasismat_eval <- Basis_eval$D1quadbasismat
    D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
    
    ## Evaluate the curves ##
    ww <- basismat%*%c[i,]
    dw <- D1quadbasismat_eval%*%c[i,]
    w <- D0quadbasismat_eval%*%c[i,]
    tau_index = which(quadpts_eval == parameters[i,3])[1]
    w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
    
    
    ##sample lambda using Gibbs
    lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, w, dw, w_delay, quadwts_eval, alpha)
    
    ##sample sigma
    sigma[i] <- GB_sigma(y, ww, alpha)
    
    ##sample c
    c_llOld <- c_ll(ww, y, sigma[i])
    c_priorOld <- c_logprior(lambda[i], parameters[i,], dw, w_delay, quadwts_eval)
    MHstep_c <- MH_c(c_priorOld, c_llOld, lambda[i], sigma[i], parameters[i,], c[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y, alpha)
    c[i,] <- MHstep_c$c
    
    ww <- basismat%*%c[i,]
    dw <- D1quadbasismat_eval%*%c[i,]
    w <- D0quadbasismat_eval%*%c[i,]
    tau_index = which(quadpts_eval == parameters[i,3])[1]
    w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
    
    
    ##sample parameters
    parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i], c[i,], dw, w, w_delay, quadwts_eval, alpha)

    
  }
  proc.time() - ptm
  
  
  parametersStore[[iter]] <- parameters
  sigmaStore[[iter]] <- sigma
  cStore[[iter]] <- c
  lambdaStore[[iter]] <- lambda
  
  
  particles <- cbind(c, parameters, sigma, lambda)
  
  print(apply(parametersStore[[iter]], 2, mean))
  print(mean(sigmaStore[[iter]]))
}


