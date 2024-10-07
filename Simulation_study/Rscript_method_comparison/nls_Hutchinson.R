library(parallel)
library(doParallel)
library(deSolve)
library(foreach)
library(dplyr)
library(Metrics)
# uncomment one of them 
spacing = 2
# spacing = 1
#  spacing = .5
# spacing = .25

dir.create("optimout/", showWarnings = FALSE)

tvec = seq(0, 30, spacing)

num_for_cores = round(detectCores()/2, digits = 0)
nsim = 300
cl = makeCluster(num_for_cores)
registerDoParallel(cl)
sigma = 0.1
tlen = length(tvec)
LVdede = function(t, y, parms) {
  r <- parms[1]
  K <- parms[2]
  tau <- parms[3]
  x_history <- parms[4]

  if (t > tau)
    lag1 = lagvalue(t - tau)
  else
    lag1 = x_history
  dx = r  * (1 - exp(lag1) / (1000 * K))
  list(c(dx))
}

out = data.frame(dede(func = LVdede, y = log(3500), times = tvec, parms = c(0.8, 2, 3, log(3500))))

objective_function <- function(parms, tvec, observed_y) {

  predicted_y = data.frame(dede(func = LVdede, y = parms[4], tvec, parms = parms))$X1
  sum((observed_y - predicted_y)^2)
}

foreach (
  i = 1:nsim,
  .errorhandling = "pass",
  .packages = c("deSolve")
) %dopar%
  {
    set.seed(i)
    if (file.exists(paste0("optimout/", i, "_", spacing, ".rds")))
      next
    
    outseed <- out
    outseed$y = outseed$X1 + rnorm(tlen, 0, sigma)
    y = outseed[, c(1, 3)]
    
    fitlist <- list()
    for (att in 1:100) {
      initial_params = c(runif(1, 1e-4, 5),runif(1, 1e-4, 5), runif(1, 1e-4, 5), y$y[1])
      fitlist[[att]] = try(optim(
        par = initial_params,
        fn = objective_function,
        tvec = tvec,
        method = "L-BFGS-B",
        lower = c(1e-4,1e-4,1e-4,5),
        upper = c(5,5,5,10),
        observed_y = y$y,
        hessian = TRUE
      ))
      
    }
    
    bestfit <- fitlist[[which.min(lapply(fitlist, function(x) ifelse(length(x) >1, x$value, Inf)))]]
    saveRDS(bestfit, file = paste0("optimout/", i, "_", spacing, ".rds"))
  }

allresult <- list()
for (i in 1:nsim) {
  allresult[[i]] <- readRDS(paste0("optimout/", i, "_", spacing, ".rds"))
}

truth = out$X1
compute_df = function(result_file) {
  df_para = data.frame(
    parameter = rep(c("r", "K", "tau", "Init"), nsim),
    est = rep(NA, nsim * 4)
  )
  df_rmse_traj = data.frame(iter = seq(1, nsim, 1),
                            rmse.traj = rep(NA, nsim))

  for (i in 1:nsim) {
    result_current = allresult[[i]]
    traj_current = rep(NA, tlen)
    df_para$est[(4 * i - 3):(4 * i )] = result_current$par
    traj_recon = dede(func = LVdede, y = df_para$est[(4 * i)], parms = c( df_para$est[(4 * i - 3):(4 * i )]),tvec)[,2]
    df_rmse_traj$rmse.traj[i] = rmse(exp(truth), exp(traj_recon))  
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
      est.var = var(est)
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
r = .8
K = 2
tau = 3
x_history = log(3500)

bias = as.numeric(unlist(tib$tib_para[, 2])) -
  c(x_history, K, r, tau)
rmse_para = sqrt(bias ^ 2 + as.numeric(unlist(
  tib$tib_para[, 4]
)))
tib_para = cbind(tib$tib_para, bias = bias, rmse = rmse_para)
knitr::kable(tib_para)
tib$rmse.traj.mean

