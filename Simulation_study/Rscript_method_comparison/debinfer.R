library(deBInfer)
library(dplyr)
library(tidyr)
library(deSolve)
library(coda)
library(Metrics)
library(foreach)
library(doParallel)
library(parallel)
library(mvtnorm)
## uncomment one of them 
# spacing = 2
# spacing = 1
# spacing = 0.5
spacing = 0.25
sigma = 0.1


tvec = seq(0, 30, spacing)

nsim = 300

DE_model = function(t, y, parms) {
  
  with(as.list(c(y, parms)), {
    if (t < tau) {
      if(t == 0) lag1 <- y
      else lag1 <- lagvalue(0)
    }
    else lag1 = lagvalue(t - tau)
    
    dy = r*(1-exp(lag1)/(1000*K))
    return(list(dy = dy))
  })
}


# the observation model: N_obs ~ normal(mu_N, var_N)
log_obs_model = function(data, sim.data, samp){
  epsilon <-1e-6
  llik.N = sum(dnorm(data$N_noisy, mean = sim.data[,"N"] + epsilon,
                     sd = samp[["sd.N"]], log = TRUE))
  return(llik.N)
}



parms = c(r = .8, K = 2, tau = 3)
N_obs = data.frame(dede(func = DE_model, parms=parms, y = log(3500), tvec))
parms["sd.N"] = sigma
colnames(N_obs) = c("time", "N")

dosim <- function(i){
  set.seed(i)
  N_obs$N_noisy = rnorm(nrow(N_obs), N_obs$N, parms["sd.N"])
  
  r = debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                   value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 5),
                   prop.var =  0.02, samp.type = "rw-ref")
  
  K = debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                   value = runif(1, 1, 2.5), prior = "unif", hypers = list(min = 0, max =5),
                   prop.var = 0.01, samp.type = "rw-ref")
  
  tau = debinfer_par(name = "tau", var.type = "de", fixed = FALSE,
                     value = runif(1, 2, 3.5), prior = "unif", hypers = list(min = 0, max = 5),
                     prop.var = 0.015, samp.type = "rw-ref")
  
  
  sd.N = debinfer_par(name = "sd.N", var.type = "obs", fixed = FALSE,
                      value = runif(1, 0, 0.2), prior = "unif", hypers = list(min = 0, max = 0.2),
                      prop.var = 0.01, samp.type = "rw-ref")
  
  N = debinfer_par(name = "N", var.type = "init", fixed = FALSE,
                   value = N_obs$N_noisy[1], prior = "unif", hypers = list(min = 5, max = 10),
                   prop.var = 0.1, samp.type = "rw-ref")
  
  mcmc.pars = setup_debinfer(r, K, tau, sd.N, N)
  
  
  mcmc_samples = de_mcmc(N = 40000, data = N_obs, de.model = DE_model,
                         obs.model = log_obs_model, all.params = mcmc.pars,
                         Tmax = max(N_obs$time), data.times = N_obs$time, cnt = 500,
                         plot = FALSE, verbose.mcmc = TRUE, solver = "dede")
  
  return(mcmc_samples)
  
}




options(mc.cores = round(detectCores()/2, digits = 0))
results <- mclapply(1:nsim, dosim)

save(results, file = paste0("debinfer_by_", spacing, ".RData"))

rmse_vec = rep(NA, nsim)
tau_vec = rep(NA, nsim)
tau_ci_width = rep(NA, nsim)
r_vec = rep(NA, nsim)
r_ci_width = rep(NA, nsim)
k_vec = rep(NA, nsim)
k_ci_width = rep(NA, nsim)
sigma_vec = rep(NA, nsim)
sigma_ci_width = rep(NA, nsim)
N_vec = rep(NA, nsim)
N_ci_width = rep(NA, nsim)
for(i in 1:nsim){
  result = results[[i]]
  sigma_vec[i] = summary(result$samples, autoburnin = TRUE)$statistics[4,1]
  sigma_ci_width[i] = summary(result$samples, autoburnin = TRUE)$quantile[4,5]- summary(result$samples, autoburnin = TRUE)$quantile[4,1]
  tau_vec[i] = summary(result$samples, autoburnin = TRUE)$statistics[3,1]
  tau_ci_width[i] =  summary(result$samples, autoburnin = TRUE)$quantile[3,5]- summary(result$samples, autoburnin = TRUE)$quantile[3,1]
  r_vec[i] = summary(result$samples, autoburnin = TRUE)$statistics[1,1]
  r_ci_width[i] = summary(result$samples, autoburnin = TRUE)$quantile[1,5]- summary(result$samples, autoburnin = TRUE)$quantile[1,1]
  k_vec[i] = summary(result$samples)$statistics[2,1]
  k_ci_width[i] =summary(result$samples, autoburnin = TRUE)$quantile[2,5]- summary(result$samples, autoburnin = TRUE)$quantile[2,1]
  N_vec[i] = summary(result$samples, autoburnin = TRUE)$statistics[5,1]
  N_ci_width[i] =summary(result$samples, autoburnin = TRUE)$quantile[5,5]- summary(result$samples, autoburnin = TRUE)$quantile[5,1]
  parms = c(r = r_vec[i], K = k_vec[i],tau = tau_vec[i])
  rmse_vec[i] = rmse(exp(dede(func = DE_model, parms= parms,
                              y = N_vec[i], tvec)[,2]), exp(N_obs$N))
  print(i)
}

df = data.frame(tau = tau_vec,
                r = r_vec,
                k = k_vec,
                sigma = sigma_vec,
                tau_ci_width = tau_ci_width,
                r_ci_width = r_ci_width,
                k_ci_width = k_ci_width,
                sigma_ci_width = sigma_ci_width,
                N = N_vec,
                N_ci_width = N_ci_width)

df_long_para = data.frame(Parameter = rep(NA, 5 * nsim),
                          Estimate = rep(NA, 5 * nsim),
                          CI_width = rep(NA, 5 * nsim))
df_long_para$Parameter = c(rep("tau",nsim),
                           rep("r",nsim),
                           rep("k",nsim),
                           rep("sigma",nsim),
                           rep("N",nsim))
df_long_para$Estimate[1: nsim] = df$tau
df_long_para$Estimate[(nsim + 1): (2*nsim)] = df$r
df_long_para$Estimate[(2*nsim + 1): (3*nsim)] = df$k
df_long_para$Estimate[(3*nsim + 1): (4*nsim)] = df$sigma
df_long_para$Estimate[(4*nsim + 1): (5*nsim)] = df$N


df_long_para$CI_width[1: nsim] = df$tau_ci_width
df_long_para$CI_width[(nsim + 1): (2*nsim)] = df$r_ci_width
df_long_para$CI_width[(2*nsim + 1): (3*nsim)] = df$k_ci_width
df_long_para$CI_width[(3*nsim + 1): (4*nsim)] = df$sigma_ci_width
df_long_para$CI_width[(4*nsim + 1): (5*nsim)] = df$N_ci_width


tib = df_long_para%>%group_by(Parameter)%>%summarise(mean = mean(Estimate),
                                                     sd= sd(Estimate),
                                                     var = var(Estimate),
                                                     mean.ci.width = mean(CI_width))

tib$bias = tib$mean - c(log(3500),2, .8, .1, 3)

tib$RMSE = sqrt((tib$bias)^2+tib$var)
knitr::kable(tib, "markdown")
mean(rmse_vec)
