library(deSolve)



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
            8 * 10 ^(-2) ## L_e

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
  
  


if (t > tau_M) Lag_M <- lagvalue(t - tau_M)[3] else Lag_M <-all_param$x0[3] 


if (t > tau_B) Lag_B <- lagvalue(t- tau_B)[1] else Lag_B <- all_param$x0[1]


if (t > tau_B + tau_P) Lag_P_add_B <- lagvalue(t-(tau_B + tau_P))[1] else Lag_P_add_B <- all_param$x0[1]   
  
  
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
                          


















