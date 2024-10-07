library(magidde)
library(parallel)
library(doParallel)
library(rstan)
library(deSolve)
library(foreach)
library(Metrics)
library(dplyr)

num_for_cores = round(detectCores() / 2, digits = 0)
nsim = 100
cl = makeCluster(num_for_cores)
registerDoParallel(cl)
sigma = 0.1
spacing = 2
## uncomment one of them to allow changes in discretization level
nlevel = 0
# nlevel = 1
# nlevel = 2
# nlevel = 3



tvec_obs = seq(0, 30, spacing)
tlen_obs = length(tvec_obs)


LVdede = function(t, y, parms){
  r <- parms[1]
  K <- parms[2]
  tau <- parms[3]
  x_history <- parms[4]
  if (t > tau) lag1 = lagvalue(t - tau) else lag1 = x_history
  dx = r  * (1 - exp(lag1)/(1000 * K))
  list(c(dx))
}


out = data.frame(dede(func = LVdede,  y = log(3500), times = tvec_obs, 
                      parms = c(0.8, 2, 3, log(3500))))


discretization_set = setDiscretization(out, nlevel)
tvec_discretization = discretization_set$time
tlen_discretization = length(discretization_set$time)
truth_discretization = data.frame(dede(func = LVdede,  y = log(3500), times = tvec_discretization, 
                      parms = c(0.8, 2, 3, log(3500))))
foo  = outer(tvec_discretization, t(tvec_discretization), '-')[, 1,]
distance_matrix = abs(foo)
sign_matrix = -sign(foo)


out_est_lag = data.frame(dede(func = LVdede,  y = log(3500),
                      parms = c(0.8, 2, 3, log(3500)), 
                      times = seq(0, 30, 0.5^(nlevel))))



model = stan_model("prey_dde_fully_bayes.stan")
allresult = foreach (
  i = 1:nsim,
  .errorhandling = "pass",
  .packages = c("magidde", "rstan")
) %dopar%
  {
    set.seed(i)
    
    yobs = out$X1 + rnorm(tlen_obs, 0, sigma)
    phi_estimate = gpsmoothing(yobs,
                               tvec_obs,
                               sigma = sigma,
                               kerneltype = "matern")
    cov = calCov(
      phi = phi_estimate$phi,
      distance_matrix,
      sign_matrix,
      bandsize = 20,
      kerneltype = "matern",
      noiseInjection = 1e-6
    )
    
    cov_c_prime = cov$Cprime
    cov_c = cov$C
    cov_c_double_prime = cov$Cdoubleprime
    
    
    if (nlevel == 0) {
      x_init = c(rep(log(3500), floor(3 / spacing) + 1), 
                 out_est_lag$X1[seq(2,length(out_est_lag$time) - floor(3 / spacing) - 1, spacing)], truth_discretization$X1)
    }else{
      x_init = c(rep(log(3500), floor(3 / (0.5^(nlevel-1))) + 1), 
                 truth_discretization$X1[seq(from = 2, by = 1, 
                                        length.out = tlen_discretization - floor(3 / (0.5^(nlevel-1)))-1)], 
                 truth_discretization$X1)
    }
    
    initf1 = function() {
      list(
        tau = 3,
        r = .8,
        k = 2,
        x_t_inc_lag = x_init
      )
    }
    
    stan_used_data = list(
      tlen_obs = tlen_obs,
      tlen_discretization = tlen_discretization,
      tvec_obs = tvec_obs,
      tvec_discretization = tvec_discretization,
      yobs = yobs,
      sigma = sigma,
      num_adjacent_pts = 2^(nlevel),
      noise_injection = 1e-6,
      
      cov_c_prime = cov_c_prime,
      cov_c = cov_c,
      cov_c_double_prime = cov_c_double_prime,
      
      phi_1 = phi_estimate$phi[1],
      phi_2 = phi_estimate$phi[2]
    )
    stan_result = sampling(
      model,
      data = stan_used_data,
      iter = 10000,
      seed = i,
      chains = 1,
      init = initf1
    )
    stan_result
  }
stopCluster(cl)
save.image(paste0("fully_level", nlevel, ".RData"))
load(paste0("fully_level", nlevel, ".RData"))

compute_infer_rmse = function(stan_fit){
  traj = rep(NA, tlen_obs)
  for(i in 1: tlen_obs){
    re = summary(stan_fit, pars = paste0('x_t_inc_lag[',2^(nlevel) * i - (2^(nlevel)-1)+tlen_discretization,']'))
    traj[i] = re$summary[1]
  }
  infer_rmse = rmse(exp(out$X1), exp(traj))
  return(list(infer_rmse = infer_rmse))
}

df_generator = function(result_file){
  df_para= data.frame(parameter = rep(c("r", "k", "tau", "N"),nsim),
                      est = rep(NA, 4 * nsim),
                      ci.width =  rep(NA, 4 * nsim),
                      sd = rep(NA, 4 * nsim))
  df_rmse_traj = data.frame(iter = seq(1, nsim,1), 
                            infer_rmse = rep(NA, nsim),
                            traj_rmse = rep(NA, nsim))
  for(i in 1:nsim){
    stan_fit = result_file[[i]]
    summarized = summary(stan_fit, pars = c('r', 'k','tau', paste0('x_t_inc_lag[',1+tlen_discretization,']')))$summary
    df_para$est[(4 * i -3 ): (4 * i)] = as.numeric(summarized[,1])
    df_para$ci.width[(4 * i -3 ): (4 * i)] = as.numeric(summarized[,8] - summarized[,4])
    df_para$sd[(4 * i -3 ): (4 * i)] = as.numeric(summarized[,3])
    df_rmse_traj$infer_rmse[i] = compute_infer_rmse(stan_fit)
    df_rmse_traj$traj_rmse[i] = rmse(exp(out$X1), 
                                     exp(data.frame(dede(func = LVdede,  
                                                         y = df_para$est[(4 * i)],
                                                         times = tvec_obs, 
                                                         parms = df_para$est[(4*i-3) : (4*i)])))$X1)
    
  }
  return(list(df_rmse_traj = df_rmse_traj,
              df_para = df_para))
}

compute_statistic = function(df_para, df_traj){
  df_para = na.omit(df_para)
  df_traj = na.omit(df_traj)
  tib = df_para %>%
    group_by(parameter) %>%
    summarise(est.mean = mean(est),
              est.sd = sd(est),
              est.var = var(est),
              ciwidth.mean = mean(ci.width),
    )
  infer_rmse_mean = mean(na.omit(unlist(df_traj$infer_rmse)))
  traj_rmse_mean = mean(na.omit(unlist(df_traj$traj_rmse)))
  return(list(tib = tib,
              infer_rmse_mean = infer_rmse_mean,
              traj_rmse_mean = traj_rmse_mean))
  
}

df_all = df_generator(allresult)
df_para = df_all$df_para
df_traj = df_all$df_rmse_traj

tib = compute_statistic(df_para, df_traj)
bias = as.numeric(unlist(compute_statistic(df_para, df_traj)$tib[,2]))-c(log(3500), 2, .8, 3)
rmse_para = sqrt(bias^2 + as.numeric(unlist(compute_statistic(df_para, df_traj)$tib[,4])))
para_tib = cbind(tib$tib, bias = bias, rmse = rmse_para)
print(tib$infer_rmse_mean)
print(tib$traj_rmse_mean)
knitr::kable(para_tib)

