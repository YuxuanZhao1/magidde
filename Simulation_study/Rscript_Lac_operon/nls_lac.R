library(deSolve)
library(dplyr)
library(Metrics)
library(doParallel)
library(parallel)
library(foreach)
dir.create("optimoutlac/", showWarnings = FALSE)
sigma = c(0.00003, 0.00001, 0.02, 0.01, 0.0005)
all_param = list(
  time =  c(seq(0, 2, 0.25), seq(3, 10, 1), seq(12, 20, 2), 25),
  theta = c(
            0.411, ## gamma_M
            0.52, ## gamma_A
            7200, ## K
            1.81 , ## K_L1
            1.95 , ## K_A
            
            0, ## gamma_L
            2880, ## alpha_L
            0.26, ## K_Le
            8.33 * 10 ^(-4), ## gamma_B
            7.25 * 10 ^(-7), ## gamma_0
            
            9.97 * 10 ^(-4), ## alpha_M
            1.76 * 10 ^(4), ## alpha_A
            1.66 * 10 ^(-2), ## alpha_B 
            2.15 * 10 ^(4), ## beta_A
            0.97, ## K_L  
            
            0.65, ## gamma_P
            10, ## alpha_P
            2.65 * 10 ^(3), ## beta_L1
            1.76 * 10 ^(4), ## beta_L2
            2.52 * 10 ^(4) ,  ## K_1 2.52 * 10 ^(-8)
            
            2,## n
            2, ## tau_B
            0.1, ## tau_M 
            0.83, ## tau_P
            2.26 * 10 ^(-2), # mu_bar
            8 * 10 ^(-2), ## L_e
            3.8 * 10 ^(-2), ## A_0
             6.26 * 10 ^(-4) ## M_0

  ),
  x0 = c(
         6.26 * 10 ^(-4), ## M_0
         0, ## B_0
         3.8 * 10 ^(-2), ## A_0
         3.72 * 10 ^(-1), ## L_0
         1.49 * 10^(-2)  ## P_0
  ))









cellmodelODE = function(t, y, theta){
  M = y[1]
  B = y[2]
  A = y[3]
  L = y[4]
  P = y[5]


  gamma_M = theta[1]
  gamma_A = theta[2]
  K = theta[3]
  K_L1 = theta[4]  
  K_A = theta[5] 
  gamma_L = theta[6] 
  alpha_L = theta[7] 
  K_Le = theta[8] 
  gamma_B = theta[9] 
  gamma_0 = theta[10] 
  alpha_M = theta[11] 
  alpha_A = theta[12] 
  alpha_B = theta[13] 
  beta_A = theta[14] 
  K_L = theta[15] 
  gamma_P = theta[16] 
  alpha_P = theta[17] 
  beta_L1 = theta[18] 
  beta_L2 = theta[19] 
  K_1 = theta[20] 
  n = theta[21] 
  tau_B = theta[22]
  tau_M = theta[23]
  tau_P = theta[24]
  mu_bar = theta[25]
  L_e = theta[26]
  A_init = theta[27]
  M_init = theta[28]

if (t > tau_M) Lag_M <- lagvalue(t - tau_M)[3] else Lag_M <- A_init 


if (t > tau_B) Lag_B <- lagvalue(t- tau_B)[1] else Lag_B <- M_init


if (t > tau_B + tau_P) Lag_P_add_B <- lagvalue(t-(tau_B + tau_P))[1] else Lag_P_add_B <- M_init   
  
  
  dM = alpha_M * (1+ K_1 * (exp(-mu_bar * tau_M) * Lag_M)^n)/
                (K+ K_1 * (exp(-mu_bar * tau_M) * Lag_M)^n) +
    gamma_0 - (gamma_M + mu_bar) * M 


  dB = alpha_B * exp(-mu_bar * tau_B) *Lag_B -(gamma_B + mu_bar) * B 
  
  dA = alpha_A * B * L/ (K_L + L) - beta_A * B * A / (K_A + A) - 
    (gamma_A +mu_bar) * A
  
  dL = alpha_L * P * L_e/(K_Le + L_e) -beta_L1 * P * L/(K_L1 + L) - 
    beta_L2 * B * L/(K_L + L) - (gamma_L + mu_bar) * L 
  
  dP = alpha_P * exp(-mu_bar * (tau_B + tau_P)) * Lag_P_add_B - (gamma_P+mu_bar) * P 
  

    
  return(list(c(dM, dB, dA, dL, dP)))
}



objective_function <- function(parms, tvec, observed_y) {
  tau_B <- parms[1]
  tau_M <- parms[2]
  tau_P <- parms[3]
  gamma_A <- parms[4]
  alpha_M <-  parms[5]
  alpha_B <- parms[6]
  alpha_P <- parms[7]
  mu_bar <- parms[8]
  M_init <- parms[9]
  B_init <- parms[10]
  A_init <- parms[11]
  L_init <- parms[12]
  P_init <- parms[13]
  


  predicted_y = data.frame(dede(func = cellmodelODE, y = c(M_init, B_init, A_init, L_init, P_init), times = tvec, parms = c(theta = c(
            0.411, ## gamma_M
            gamma_A, ## gamma_A
            7200, ## K
            1.81 , ## K_L1
            1.95 , ## K_A
            
            0, ## gamma_L
            2880, ## alpha_L
            0.26, ## K_Le
            8.33 * 10 ^(-4), ## gamma_B
            7.25 * 10 ^(-7), ## gamma_0
            
            alpha_M, ## alpha_M
            1.76 * 10 ^(4), ## alpha_A
            alpha_B, ## alpha_B 
            2.15 * 10 ^(4), ## beta_A
            0.97, ## K_L  
            
            0.65, ## gamma_P
            alpha_P, ## alpha_P
            2.65 * 10 ^(3), ## beta_L1
            1.76 * 10 ^(4), ## beta_L2
            2.52 * 10 ^(4) ,  ## K_1 2.52 * 10 ^(-8)
            
            2,## n
            tau_B, ## tau_B
            tau_M, ## tau_M 
            tau_P, ## tau_P
            mu_bar, # mu_bar
            8 * 10 ^(-2), ## L_e
            A_init,
            M_init


  )), atol = 1e-5, rtol = 1e-5))
  #  rtol and atol: relative and absolute tolerances that define error control. 
  # Default is 1e-6. This aviods from be stuck in numerically solving the lac operon system. 
  1/(sigma[1])^2 * sum((observed_y$noise_X1 - predicted_y$X1)^2) +   1/(sigma[2])^2 *sum((observed_y$noise_X2 - predicted_y$X2)^2) + 1/(sigma[3])^2 * sum((observed_y$noise_X3 - predicted_y$X3)^2) +  1/(sigma[4])^2 *  sum((observed_y$noise_X4 - predicted_y$X4)^2) +   1/(sigma[5])^2 * sum((observed_y$noise_X5 - predicted_y$X5)^2)
}



generated_data_dde = data.frame(dede(
  func = cellmodelODE,
  y = all_param$x0,
  times = all_param$time,
  parms = all_param$theta
))



num_for_cores = round(detectCores() /2, digits = 0)
nsim = 100
cl = makeCluster(num_for_cores)
registerDoParallel(cl)
allresult =
  foreach (
    i = 1:nsim,
    .errorhandling = "pass",
    .packages = c("deSolve")
  ) %dopar%
  {
    set.seed(i)
    if (file.exists(paste0("optimoutlac/", i,  ".rds")))
      next
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
    fitlist <- list()
    for (att in 1 : 100) {
  initial_params = c(runif(1, 0, 5),runif(5, 0, 1), runif(1, 8, 12), runif(1, 0, 1), as.numeric(y[1,2:6]))
    fitlist[[att]] = optim(
      par = initial_params,
      fn = objective_function,
      tvec = all_param$time,
      method = "L-BFGS-B",
        lower = c(rep(0, 6), 8, rep(0, 6)),
      upper = c(5, rep(1, 5), 12, rep(1, 6)),
      observed_y = y,
      # ndeps: A vector of step sizes for the finite-difference approximation to the gradient, default is 1e-3. 
      # Here we need a small ndeps due to small parameter values. 
      control = list(ndeps = c(rep(1e-5, 6), 1e-3, rep(1e-7, 6)))
    )}
    

    bestfit <- fitlist[[which.min(lapply(fitlist, function(x) ifelse(length(x) >1, x$value, Inf)))]]
    saveRDS(bestfit, file = paste0("optimoutlac/", i, ".rds"))
  }

allresult <- list()
for (i in 1:nsim) {
  allresult[[i]] <- readRDS(paste0("optimoutlac/", i, ".rds"))
}

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
  est = rep(NA, nsim * 13)
)

df_traj = data.frame(component = rep(c("x1","x2","x3","x4","x5"), nsim),
                     rmse = rep(NA, 5 * nsim))

for (i in 1:nsim) {
  result = allresult[[i]]
  coef = result$par
  df_para$est[(13 * i - 12):(13 * i)] = coef
  predicted_y = data.frame(dede(func = cellmodelODE, y = c(coef[9], coef[10], coef[11], coef[12], coef[13]), times = all_param$time, parms = c(theta = c(
            0.411, ## gamma_M
            coef[4], ## gamma_A
            7200, ## K
            1.81 , ## K_L1
            1.95 , ## K_A
            
            0, ## gamma_L
            2880, ## alpha_L
            0.26, ## K_Le
            8.33 * 10 ^(-4), ## gamma_B
            7.25 * 10 ^(-7), ## gamma_0
            
            coef[5], ## alpha_M
            1.76 * 10 ^(4), ## alpha_A
            coef[6], ## alpha_B 
            2.15 * 10 ^(4), ## beta_A
            0.97, ## K_L  
            
            0.65, ## gamma_P
            coef[7], ## alpha_P
            2.65 * 10 ^(3), ## beta_L1
            1.76 * 10 ^(4), ## beta_L2
            2.52 * 10 ^(4) ,  ## K_1 2.52 * 10 ^(-8)
            
            2,## n
            coef[1], ## tau_B
            coef[2], ## tau_M 
            coef[3], ## tau_P
            coef[8], # mu_bar
            8 * 10 ^(-2), ## L_e
            coef[11],
            coef[9]))))
  df_traj$rmse[5 * i -4] = rmse(predicted_y$X1, generated_data_dde$X1)
  df_traj$rmse[5 * i -3] = rmse(predicted_y$X2, generated_data_dde$X2)
  df_traj$rmse[5 * i -2] = rmse(predicted_y$X3, generated_data_dde$X3)
  df_traj$rmse[5 * i -1] = rmse(predicted_y$X4, generated_data_dde$X4)
  df_traj$rmse[5 * i] = rmse(predicted_y$X5, generated_data_dde$X5)
  print(i)
}

tib_para = df_para %>% group_by(param) %>% summarise(
  est.mean = mean(est),
  est.sd = sd(est),
  est.var = var(est)
)

bias = as.numeric(unlist(tib_para[, 2] - c(all_param$theta[c(22, 23, 24, 2, 11, 13, 17, 25)], all_param$x0)))
rmse = as.numeric(unlist(sqrt(bias ^ 2 + tib_para[, 4])))
cbind(tib_para, bias = bias, rmse = rmse)
knitr::kable(cbind(tib_para, bias = bias, rmse = rmse))


tib_traj = df_traj %>% group_by(component) %>% summarise(
  est.mean = mean(rmse),
  est.sd = sd(rmse)
)

knitr::kable(tib_traj)
