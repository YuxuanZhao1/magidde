library(magidde)
library(parallel)
library(doParallel)
library(deSolve)
library(foreach)
library(dplyr)
library(Metrics)
# uncomment one of them 
spacing = 2
# spacing = 1
# spacing = .5
# spacing = .25


tvec = seq(0, 30, spacing)

num_for_cores = round(detectCores()/2, digits = 0)
nsim = 300
cl = makeCluster(num_for_cores)
registerDoParallel(cl)
sigma = 0.1
spacing = 2
## uncomment one of them to allow changes in discretization level
nlevel = 0
# nlevel = 1
# nlevel = 2
# nlevel = 3



tvec = seq(0, 30, spacing)
tlen = length(tvec)


LVdede = function(t, y, parms){
  r <- parms[1]
  K <- parms[2]
  tau <- parms[3]
  x_history <- parms[4]
  if (t > tau) lag1 = lagvalue(t - tau) else lag1 = x_history
  dx = r  * (1 - exp(lag1)/(1000 * K))
  list(c(dx))
}


out = data.frame(dede(func = LVdede,  y = log(3500), times = tvec, 
                      parms = c(0.8, 2, 3, log(3500))))



allresult =
  foreach (
    i = 1:nsim,
    .errorhandling = "pass",
    .packages = c("magidde")
  ) %dopar%
  {
    set.seed(i)
    out$y = out$X1 + rnorm(tlen, 0, sigma)
    y = out[, c(1, 3)]
    if (nrow(y) <= 61) {
      # set discretization level by 0.5 for |gamma| = 16, 31, 61
      yIn = setDiscretization(y, by = 0.5)
    } else {
      # discretization level is same as observations for |gamma| = 121
      yIn = y
    }
    lvmodel = list(
      name = "lvD",
      thetaLowerBound = c(0, 0, 0),
      thetaUpperBound = c(Inf, Inf, Inf)
    )
    res1 =
      MagiSolver(yIn,
                 lvmodel,
                 control = list(nstepsHmc = 20,
                                niterHmc = 40000, kerneltype = "matern"))
    res1
  }
stopCluster(cl)
save.image(paste0("magi_by", spacing, ".RData"))

load(file = paste0("magi_by", spacing, ".RData"))

compute_df = function(result_file) {
  df_para = data.frame(
    parameter = rep(c("r", "K", "tau", "sigma", "Init"), nsim),
    est = rep(NA, nsim * 5),
    ci.width =  rep(NA, nsim * 5)
  )
  df_rmse_traj = data.frame(iter = seq(1, nsim, 1),
                            rmse.traj = rep(NA, nsim),
                            rmse.infer = rep(NA, nsim))
  for (i in 1:nsim) {
    result_current = allresult[[i]]
    traj_current = rep(NA, tlen)
    df_para$est[(5 * i - 4):(5 * i - 1)] = as.numeric(summary(result_current, sigma = TRUE)[1, ])
    df_para$ci.width[(5 * i - 4):(5 * i - 1)] = as.numeric(summary(result_current, sigma = TRUE)[3, ] - summary(result_current, sigma = TRUE)[2, ])
    df_para$est[(5 * i)] = colMeans(result_current$xsampled[, , 1])[1]
    df_para$ci.width[5 * i] =  as.numeric(
      quantile(result_current$xsampled[, , 1][, 1], 0.975) - quantile(result_current$xsampled[, , 1][, 1], 0.025)
    )

    traj_recon = dede(func = LVdede, times = tvec, 
                      parms = c(df_para$est[c(5 * i - 4,
                                              5 * i - 3,
                                              5 * i - 2,
                                              5 * i)]),
                      y =  df_para$est[5 * i])[,2]
    df_rmse_traj$rmse.traj[i] = rmse(exp(out$X1), exp(traj_recon))  
  
  }
  return(list(df_para = df_para,
              df_rmse_traj = df_rmse_traj))
}

compute_statistic = function(df_para, df_traj) {
  tib_para = df_para %>%
    group_by(parameter) %>%
    summarise(
      est.mean = mean(est),
      est.sd = sd(est),
      est.var = var(est),
      ciwidth.mean = mean(ci.width),
    )

  rmse.traj.mean = mean(df_traj$rmse.traj)
  rmse.traj.sd = sd(df_traj$rmse.traj)
  return(list(
    tib_para = tib_para,
    rmse.traj.mean = rmse.traj.mean,
    rmse.traj.sd = rmse.traj.sd
  ))
  
}

df_all = compute_df(allresult)
df_para = df_all$df_para
df_traj = df_all$df_rmse_traj
tib = compute_statistic(df_para, df_traj)

bias = as.numeric(unlist(tib$tib_para[, 2])) -
  c(log(3500), 2, .8, sigma, 3)
rmse_para = sqrt(bias ^ 2 + as.numeric(unlist(
  tib$tib_para[, 4]
)))
tib_para = cbind(tib$tib_para, bias = bias, rmse = rmse_para)
knitr::kable(tib_para)
tib$rmse.traj.mean

