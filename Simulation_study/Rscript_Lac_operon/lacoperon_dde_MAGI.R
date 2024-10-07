library(magidde)
library(deSolve)
library(foreach)
library(doParallel)
library(parallel)
library(knitr)
library(dplyr)
library(Metrics)
source("lactose_operon_dde_simulate_data.R")
generated_data_dde = data.frame(dede(
  func = cellmodelODE,
  y = all_param$x0,
  times = all_param$time,
  parms = all_param$theta
))
sigma = c(0.00003, 0.00001, 0.02, 0.01, 0.0005)
num_for_cores = round(detectCores() /2, digits = 0)
nsim = 100
cl = makeCluster(num_for_cores)
registerDoParallel(cl)
allresult =
  foreach (
    i = 1:nsim,
    .errorhandling = "pass",
    .packages = c("magidde")
  ) %dopar%
  {
    set.seed(i)
    generated_data_dde$noise_X1 =  generated_data_dde$X1  + rnorm(length(generated_data_dde$X1),
                                                                  0,
                                                                  sigma[1])
    generated_data_dde$noise_X2 =  generated_data_dde$X2 + rnorm(length(generated_data_dde$X2),
                                                                 0,
                                                                 sigma[2])
    
    generated_data_dde$noise_X3 =  generated_data_dde$X3 + rnorm(length(generated_data_dde$X3),
                                                                 0,
                                                                 sigma[3])
    
    generated_data_dde$noise_X4 =  generated_data_dde$X4 + rnorm(length(generated_data_dde$X4),
                                                                 0,
                                                                 sigma[4])
    
    
    generated_data_dde$noise_X5 =  generated_data_dde$X5 + rnorm(length(generated_data_dde$X5),
                                                                 0,
                                                                 sigma[5])
    
    y = generated_data_dde[, c(1, 7:11)]
    
    lacmodel = list(
      thetaLowerBound = rep(0, 8),
      thetaUpperBound = rep(Inf, 8),
      name = "lacoperonwpars"
    )
    
    y_In = setDiscretization(y, by = .25)
    
    res_allpars = MagiSolver(y_In,
                             lacmodel,
                             control = list(
                               niterHmc = 50000,
                               sigma = sigma,
                               useFixedSigma = TRUE
                             ))
    res_allpars
  }

stopCluster(cl)
save.image(file = "lac_operon.RData")
load("lac_operon.RData")
tvec = seq(0, 25, .25)
tlen = length(tvec)

theoretical_traj = data.frame(dede(
  func = cellmodelODE,
  y = all_param$x0,
  times = all_param$time,
  parms = all_param$theta
))

df_para = data.frame(
  param = rep(
    c(
      "theta_1",
      "theta_2",
      "theta_3",
      "theta_4",
      "theta_5",
      "theta_6",
      "theta_7",
      "theta_8",
      "x1_init",
      "x2_init",
      "x3_init",
      "x4_init",
      "x5_init"
    ),
    nsim
  ),
  est = rep(NA, nsim * 13),
  ci.width = rep(NA, nsim * 13)
)

df_traj = data.frame(component = rep(c("x1","x2","x3","x4","x5"), nsim),
                     rmse = rep(NA, 5 * nsim))



for (i in 1:nsim) {
  result = allresult[[i]]
  df_para$est[(13 * i - 12):(13 * i - 5)] = as.numeric(summary(result)[1, ])
  df_para$ci.width[(13 * i - 12):(13 * i - 5)] =  as.numeric(summary(result)[3, ]) -  as.numeric(summary(result)[2, ])
  df_para$est[(13 * i - 4):(13 * i)] = as.numeric(colMeans(result$xsampled[, 1, 1:5]))
  df_para$ci.width[(13 * i - 12):(13 * i - 5)] =  as.numeric(summary(result)[3, ]) -  as.numeric(summary(result)[2, ])
  df_para$ci.width[13 * i - 4] = as.numeric(quantile(result$xsampled[, 1, 1], 0.975) -  quantile(result$xsampled[, 1, 1], 0.025))
  df_para$ci.width[13 * i - 3] = as.numeric(quantile(result$xsampled[, 1, 2], 0.975) -  quantile(result$xsampled[, 1, 2], 0.025))
  df_para$ci.width[13 * i - 2] = as.numeric(quantile(result$xsampled[, 1, 3], 0.975) -  quantile(result$xsampled[, 1, 3], 0.025))
  df_para$ci.width[13 * i - 1] = as.numeric(quantile(result$xsampled[, 1, 4], 0.975) -  quantile(result$xsampled[, 1, 4], 0.025))
  df_para$ci.width[13 * i] = as.numeric(quantile(result$xsampled[, 1, 5], 0.975) -  quantile(result$xsampled[, 1, 5], 0.025))
  df_traj$rmse[5 * i - 4] = rmse(colMeans(result$xsampled[,match(all_param$time , seq(0, 25, .25)), 1]), theoretical_traj$X1)
  df_traj$rmse[5 * i - 3] = rmse(colMeans(result$xsampled[,match(all_param$time , seq(0, 25, .25)), 2]), theoretical_traj$X2)
  df_traj$rmse[5 * i - 2] = rmse(colMeans(result$xsampled[,match(all_param$time , seq(0, 25, .25)), 3]), theoretical_traj$X3)
  df_traj$rmse[5 * i - 1] = rmse(colMeans(result$xsampled[,match(all_param$time , seq(0, 25, .25)), 4]), theoretical_traj$X4)
  df_traj$rmse[5 * i] = rmse(colMeans(result$xsampled[,match(all_param$time , seq(0, 25, .25)), 5]), theoretical_traj$X5)
  print(i)
}

tib_para = df_para %>% group_by(param) %>% summarise(
  est.mean = mean(est),
  ci.width.mean = mean(ci.width),
  est.sd = sd(est),
  est.var = var(est)
)

bias = as.numeric(unlist(tib_para[, 2] - c(all_param$theta[c(22, 23, 24, 2, 11, 13, 17, 25)], all_param$x0)))
rmse = as.numeric(unlist(sqrt(bias ^ 2 + tib_para[, 5])))
cbind(tib_para, bias = bias, rmse = rmse)
kable(cbind(tib_para, bias = bias, rmse = rmse))


tib_traj = df_traj %>% group_by(component) %>% summarise(
  est.mean = mean(rmse),
  est.sd = sd(rmse)
)

kable(tib_traj)


