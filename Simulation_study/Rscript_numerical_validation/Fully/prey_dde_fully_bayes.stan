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
  
  
  int signnum(real x) { return x < 0 ? -1 : x > 0; }
  
  matrix calCov_deriv(real[] tvec_1, real[] tvec_2, real phi_1, real phi_2, int tlen){
    
    matrix[tlen,tlen] deriv_cov_mat; 
    real dist_num;
    real dist_abs;
    for (i in 1 : tlen){
      for (j in 1 : tlen){
        if ((i == 1) || (j == 1)) {
          dist_num = tvec_1[i]-tvec_2[j];
          dist_abs = fabs(dist_num);
          deriv_cov_mat[i,j] =   -signnum(dist_num)* (phi_1 * exp((-sqrt(5)*dist_abs)/phi_2)) * (((5*dist_abs)/(3*phi_2^2)) + ((5*sqrt(5)*dist_abs^2)/(3*phi_2^3)));
        } else {
          deriv_cov_mat[i,j] = deriv_cov_mat[i-1, j-1];
        }
      }
    }
    return deriv_cov_mat;
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
    real noise_injection;
    
    matrix[tlen_discretization, tlen_discretization] cov_c_double_prime; // covariance matrix of C''
    matrix[tlen_discretization, tlen_discretization] cov_c_prime; // covariance matrix of C'
    matrix[tlen_discretization, tlen_discretization] cov_c; // covariance matrix of C
    
    real phi_1;
    real phi_2;
  }
  
  
  parameters {
    real<lower = 0> tau;
    real<lower = 0> r;
    real<lower = 0> k;
    vector[2 * tlen_discretization] x_t_inc_lag; 
  }
  
  
  
  model {
    vector[tlen_discretization] f;
    real x_real[tlen_discretization];
    real x_lag_real[tlen_discretization];
    vector[tlen_discretization] x_lag;
    real tvec_lag[tlen_discretization];
    vector[tlen_discretization] mu_x_derivative;
    matrix[tlen_discretization, tlen_discretization] sigma_x_derivative;
    
    matrix[tlen_discretization, tlen_discretization] sigma_12; 
    matrix[2 * tlen_discretization, 2 * tlen_discretization] Sigma;
    matrix[tlen_discretization, 2 * tlen_discretization] pre_mult;
    matrix[tlen_discretization, 2 * tlen_discretization] sigma_12_append; 
    
    target += uniform_lpdf(tau | 0, 5);
    target += uniform_lpdf(r |0, 5);
    target += uniform_lpdf(k |0, 5);
    
    for(i in 1:tlen_discretization){
      tvec_lag[i] = tvec_discretization[i] - tau;
    }
    
    sigma_12 = calCov(tvec_lag, tvec_discretization, phi_1, phi_2, tlen_discretization);
    Sigma = append_row(append_col(cov_c, sigma_12), append_col(sigma_12', cov_c));

    target += multi_normal_lpdf(x_t_inc_lag | rep_vector(0.0, 2* tlen_discretization), Sigma);
    
    for (i in 1:tlen_obs){
      real yi = yobs[i];
      target += normal_lpdf( x_t_inc_lag[tlen_discretization+num_adjacent_pts * i - (num_adjacent_pts-1)] | yi, sigma); 
    }
    
    for (i in 1:tlen_discretization){
      x_real[i] =  x_t_inc_lag[tlen_discretization + i];
      x_lag_real[i] =  x_t_inc_lag[i];
    }
    
    sigma_12_append = append_col(calCov_deriv(tvec_discretization, tvec_lag, phi_1, phi_2, tlen_discretization), cov_c_prime);
    
    pre_mult = mdivide_right_spd(sigma_12_append, Sigma);
    
    mu_x_derivative =  pre_mult * x_t_inc_lag;
    
    sigma_x_derivative = cov_c_double_prime - pre_mult * sigma_12_append';
    f = dde(tlen_discretization, tvec_discretization, tau, r, k, x_real,x_lag_real);
    target += multi_normal_lpdf(f|mu_x_derivative, add_diag(sigma_x_derivative, noise_injection));
    
  }