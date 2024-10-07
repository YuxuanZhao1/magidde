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
  
  real interp_1d_linear(real[] tvec, real[] x_t, real x_out) {
  int left = 1;
  int right = 1;
  real w = 1.0;

  real y_out;
   
    while (tvec[right] < x_out) {
      right = right + 1;
    }
    while (tvec[left + 1] < x_out) {
      left = left + 1;
    }
    w = (tvec[right] - x_out) / (tvec[right] - tvec[left]);
    y_out = w * x_t[left] + (1 - w) * x_t[right];
     return y_out;
  }

  
  vector dde(int tlen, real[] tvec, real tau, real r, real k, real[] x_t){
    vector[tlen] dx_t;
    vector[tlen] phase_type;
    real x_1_t;
    real tout;
    
   for (i in 1:tlen){
    phase_type[i] = phase(tvec[i], tau);
    if (phase_type[i] == 1){
      x_1_t = x_t[1];

    }else{
      tout = tvec[i] - tau;
      x_1_t = interp_1d_linear(tvec, x_t, tout);
      }
      dx_t[i] = r * ( 1 - exp(x_1_t)/(1000 * k));
   }

  return dx_t;
  }
  
  
  }
  
  
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
    real tvec_lag[tlen_discretization];
    
    target += uniform_lpdf(tau | 0, 5);
    
    target += uniform_lpdf(r | 0, 5);
    
    target += uniform_lpdf(k | 0, 5);
    
    target += multi_normal_lpdf(x_t | rep_vector(0.0,  tlen_discretization), cov_c);
    
    
    
      for (i in 1:tlen_obs){
      real yi = yobs[i];
       target += normal_lpdf( x_t[num_adjacent_pts * i - (num_adjacent_pts-1)] | yi, sigma); 
      
    }
    
    for (i in 1 : tlen_discretization){
      x_real[i] =  x_t[i];
      }
    
    f = dde(tlen_discretization, tvec_discretization, tau, r, k, x_real);
    target += multi_normal_lpdf(f| cov_m * x_t , cov_k);
  }
