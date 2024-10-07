rm(list=ls())
library("fda")
library("truncnorm")
library(numDeriv)
library(fda)
library(deSolve)
library(NLRoot)
library(Hmisc)

NITER <- 300
times <- seq(0, 30, 2)

knots = seq(min(times), max(times),by=2)
sigma1 <- 0.1
sigma1_simulator <- sigma1
seeds <- 1:NITER



alambda <- 1
blambda <- 1
sigmac <- 0.3

CESSthresholds <- 0.9
NP <- 500

filename_nu = './nu'
filename_P = './P'
filename_tau = './tau'
filename_W = './W'
filename_sigma = './sigma'
filename_MSE = './MSE'
filename_lambda = './lambda'

filename_nuCI = './nuCI'
filename_PCI = './PCI'
filename_tauCI = './tauCI'
filename_WCI = './WCI'
filename_sigmaCI = './sigmaCI'

filename_iter = './iterations'

ifAppend=FALSE;
for(iter in 1:NITER){
  if(iter == 1){
    ifAppend=FALSE
  }else{
    ifAppend=TRUE
  }
  set.seed(seeds[iter])
  source('mainDDE.R')
  mean_para <- apply(parameters, 2, mean)
  mean_r <- mean_para[1]
  mean_P <- mean_para[2]
  mean_tau <- mean_para[3]
  meansigma <- mean(sigma)
  fitted1 <- colMeans(c %*% t(basismat))
  cMean <- colMeans(c)

  quantile_para <- apply(parameters, 2, quantile, c(0.025,0.975))
  r_Credible <- as.numeric(apply(parameters, 2, quantile, c(0.025,0.975))[,1])
  P_Credible <- as.numeric(apply(parameters, 2, quantile, c(0.025,0.975))[,2])
  tau_Credible <- as.numeric(apply(parameters, 2, quantile, c(0.025,0.975))[,3])
  sigma_Credible <- as.numeric(quantile(sigma, c(0.025, 0.975)))
  W_Credible <- as.numeric(quantile((c %*% t(basismat))[,1], c(0.025, 0.975)))


  write.table(matrix(round(mean_r, digits = 5),nr=1), file = paste(filename_nu,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mean_P, digits = 5),nr=1), file = paste(filename_P,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mean_tau, digits = 5),nr=1), file = paste(filename_tau,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meansigma, digits = 5),nr=1), file = paste(filename_sigma,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(lambda, digits = 5),nr=1), file = paste(filename_lambda,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(sqrt(sum((output[,2] - exp(fitted1))^2)/length(fitted1)), digits = 5),nr=1), file = paste(filename_MSE,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(cMean[1], digits = 5),nr=1), file = paste(filename_W,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);

  write.table(matrix(round(r_Credible, digits = 5),nr=1), file = paste(filename_nuCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(P_Credible, digits = 5),nr=1), file = paste(filename_PCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(tau_Credible, digits = 5),nr=1), file = paste(filename_tauCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(sigma_Credible, digits = 5),nr=1), file = paste(filename_sigmaCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(W_Credible, digits = 5),nr=1), file = paste(filename_WCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  # write.table(matrix(round(iter, digits = 5),nr=1), file = paste(filename_iter,".txt",sep=""), append = ifAppend, row.names = FALSE,
  #             col.names = FALSE);

}



