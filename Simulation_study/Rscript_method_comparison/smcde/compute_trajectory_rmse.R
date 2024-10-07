library(Metrics)
library(deSolve)

nsim = 300
## Uncomment one of them 
spacing = 2
# spacing = 1
# spacing = .5
# spacing = .25
sigma = 0.1
tvec = seq(0, 30, spacing)

LVdede = function(t, y, parms){
  r <- parms[1]
  K <- parms[2]
  tau <- parms[3]
  x_history <- parms[4]
  if (t > tau) lag1 = lagvalue(t - tau) else lag1 = x_history
  dx = r  * (1 - exp(lag1)/(1000 * K))
  list(c(dx))
}


N_true = data.frame(dede(func = LVdede,  y = log(3500), times = tvec, parms = c(0.8, 2, 3, log(3500))))$X1

tau_vec = rep(NA, nsim)
r_vec = rep(NA, nsim)
k_vec = rep(NA, nsim)
N_vec = rep(NA, nsim)
traj_rmse = rep(NA, nsim)

tau_vec = read.table(paste0("01_", spacing, "/tau.txt"), header = FALSE)$V1
r_vec = read.table(paste0("01_", spacing, "/nu.txt"), header = FALSE)$V1
k_vec = read.table(paste0("01_", spacing, "/P.txt"), header = FALSE)$V1
N_vec= read.table(paste0("01_", spacing, "/W.txt"), header = FALSE)$V1

for(i in 1 : nsim){
  traj_current =  dede(func = LVdede, y = N_vec[i], tvec, parms = c(r_vec[i],
                                                                    k_vec[i],
                                                                    tau_vec[i],
                                                                    N_vec[i]))[,2]
  traj_rmse[i] = rmse(exp(traj_current), exp(N_true))
}

mean(traj_rmse)
write.csv(traj_rmse, paste0("01_", spacing, "trajectory_rmse.txt"))
