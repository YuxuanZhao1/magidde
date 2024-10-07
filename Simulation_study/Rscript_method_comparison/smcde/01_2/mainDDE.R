source("splineDDE.R")
source('simulatorDDE.R')

set.seed(364)
M <- nknots

r_particles = abs(rnorm(NP, 1.2, 1))
P_particles = abs(rnorm(NP, 2, 1))
tau_particles = rnorm(NP, 2, 0.6)

sigma <- runif(NP, 0.3, 5)

lambda <- (rgamma(NP, 1, 1))

c_LSE <- solve(as.matrix(basismat2)+ diag(0.000001, nrow = dim(basismat2)[1], ncol = dim(basismat2)[2] ))%*%t(basismat)%*%log(y)


#sum((output[,2]-basismat%*%c1_LSE)^2)
double_quadpts <- table(quadpts)[which(table(quadpts) == 2)]



c <- matrix(NA, nr = NP, nc = ncol(basismat))

for(i in 1:NP){
  c[i,] <- rnorm(length(c_LSE), c_LSE, 0.2)
}

resampleThreshold <- 0.5

source('functionDDE.R')
source('resample.R')
source('asmcDDE.R')
