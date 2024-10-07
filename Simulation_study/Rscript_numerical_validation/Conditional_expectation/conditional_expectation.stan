functions {
  int phase(real t, real tau){
    int i;
    
    if (t > tau ){
      i = 2;
    } 
    else{
      i = 1;
    }
    
    return i;
  }
  
  
  matrix calCov(real[] tvec_1, real[] tvec_2, real phi_1, real phi_2, int tlen){
    
    matrix[tlen,tlen] cov_mat; 
    real dist_val;
    for (i in 1 : tlen){
      for (j in 1 : tlen){
        if ((i == 1) || (j == 1)) {
          dist_val = fabs(tvec_1[i]-tvec_2[j]);
          cov_mat[i,j] = phi_1 * (1 + sqrt(5) * dist_val / phi_2 + 5 * dist_val^2/
          (3 * phi_2 ^2)) * exp(-sqrt(5) * dist_val/phi_2);
        } else {
          cov_mat[i,j] = cov_mat[i-1, j-1];
        }
      }
    }
    return cov_mat;
  }
  
  vector dde(int tlen, real[] tvec, real tau, real r, real k, real[] x_t, real[] x_t_lag){
    vector[tlen] dx_t;
    vector[tlen] phase_type;
    real x_1_t;
    
    
    for (i in 1:tlen){
      phase_type[i] = phase(tvec[i], tau);
      if (phase_type[i] == 1){
        x_1_t = x_t[1];
        
      }else{
        x_1_t = x_t_lag[i];
      }
      dx_t[i] = r * ( 1 - exp(x_1_t)/(1000 * k));
    }
    
    return dx_t;
  }}
  
  
  data {
    int tlen_obs;
    int tlen_discretization;
    real tvec_obs[tlen_obs]; 
    real tvec_discretization[tlen_discretization]; 
    real yobs[tlen_obs]; 
    real sigma;
    int num_adjacent_pts;
    
  matrix[tlen_discretization, tlen_discretization] cov_m; // covariance matrix of m
  matrix[tlen_discretization, tlen_discretization] cov_c; // covariance matrix of C
  matrix[tlen_discretization, tlen_discretization] cov_k; // covariance matrix of k
  matrix[tlen_discretization, tlen_discretization] cov_cinv;
    
    real phi_1;
    real phi_2;
    
    
  }
  
  
  parameters {
    real<lower = 0> tau;
    real<lower = 0> r;
    real<lower = 0> k;
    vector[tlen_discretization] x_t; 
  }
  
  
  
  model {
    vector[tlen_discretization] f;
    real x_real[tlen_discretization];
    real x_lag_real[tlen_discretization];
    vector[tlen_discretization] x_lag;
    real tvec_lag[tlen_discretization];
    matrix[tlen_discretization, tlen_discretization] sigma_12; 
    
    target += uniform_lpdf(tau | 0, 5);
    
    target += uniform_lpdf(r | 0, 5);
    
    target += uniform_lpdf(k | 0, 5);
    
    for(i in 1 : tlen_discretization){
      tvec_lag[i] = tvec_discretization[i] - tau;
    }
    sigma_12 = calCov(tvec_lag, tvec_discretization, phi_1, phi_2, tlen_discretization);
    
    target += multi_normal_lpdf(x_t | rep_vector(0.0,  tlen_discretization), cov_c);
    
    x_lag = (sigma_12 * cov_cinv * x_t);
    
    
      for (i in 1:tlen_obs){
      real yi = yobs[i];
       target += normal_lpdf( x_t[num_adjacent_pts * i - (num_adjacent_pts-1)] | yi, sigma); 
      
    }
    
    for (i in 1 : tlen_discretization){
      x_real[i] =  x_t[i];
      x_lag_real[i] = x_lag[i];
      }
    
    f = dde(tlen_discretization, tvec_discretization, tau, r, k, x_real, x_lag_real);
    target += multi_normal_lpdf(f | cov_m * to_vector(x_real), cov_k);
  }
