data {
  int<lower=1> N;
  int<lower=1> T;
  real std_dev[N];
  int<lower=1, upper=T> Tsubj[N];
  int<lower=0, upper=1> choice[N, T];
  real delay[N, T];
}


parameters {
  // Hyper(group)-parameters
  vector[4] mu; 
  vector<lower=0>[4] sigma;
  
  
  // individual paramereters
  vector[N] m_min_raw;
  vector[N] m_plus_raw;
  vector<lower = 0>[N] s_min_raw;
  vector<lower = 0>[N] s_plus_raw;
}

transformed parameters {
  // for Matt trick
  vector[N]  m_min;
  vector[N]  m_plus;
  vector<lower = 0>[N]  s_min;
  vector<lower=  0>[N]  s_plus;
  m_min = mu[1] + sigma[1]*m_min_raw; //   m_min_raw~normal(0,1) =>  m_min~normal(mu[1],sigma[1])
  m_plus = mu[2] + sigma[2]*m_plus_raw;  //  m_plus_raw~normal(0,1) =>  m_plus~normal(mu[2],sigma[2])
  s_min = Phi_approx(mu[3] + sigma[3]*s_min_raw); //  s_min_raw~normal(0,1) => s_min \in [0,1]
  s_plus = Phi_approx(mu[4] + sigma[4]*s_plus_raw); //  s_plus_raw~normal(0,1) => s_plus \in [0,1]
}


model {
  // Hyper-priors
  mu[1]  ~ normal(-10, 10);
  mu[2]  ~ normal(10, 10);
  mu[3]  ~ normal(0, 1);
  mu[4]  ~ normal(0, 1);
  sigma ~ student_t(4,0,1);
  
  // individual parameters for Matt Trick
  m_min_raw  ~ normal(0, 1);
  m_plus_raw  ~ normal(0, 1);  
  s_min_raw  ~ normal(0, 1);
  s_plus_raw  ~ normal(0, 1);


  for (i in 1:N) {
    real theta;
    for (t in 1:Tsubj[i]) {
      theta =   Phi_approx((m_plus[i] - delay[i,t])/(sqrt(s_plus[i]^2+std_dev[i]^2))) - Phi_approx((m_min[i] - delay[i,t])/(sqrt(s_min[i]^2+std_dev[i]^2)));
          choice[i,t] ~ bernoulli(theta);
    }
  }
}



generated quantities {
   // For group level parameters
  real mu_m_min; 
  real mu_m_plus;
  real mu_s_min;
  real mu_s_plus;
  real sigma_m_min; 
  real sigma_m_plus;
  real sigma_s_min;
  real sigma_s_plus;
  
  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N,T];

  // Set all posterior predictions to -1 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }
  mu_m_min = mu[1]; 
  mu_m_plus = mu[2];
  mu_s_min = Phi_approx(mu[3]);
  mu_s_plus = Phi_approx(mu[4]);
  sigma_m_min = sigma[1]; 
  sigma_m_plus = sigma[2];
  sigma_s_min = sigma[3];
  sigma_s_plus = sigma[4];
  
   { // local section, this saves time and space
    for (i in 1:N) {
      real theta;
      
      // Initialize values
      log_lik[i] = 0;
      for (t in 1:Tsubj[i]) {
	theta =   Phi_approx((m_plus[i] - delay[i,t])/(sqrt(s_plus[i]^2+std_dev[i]^2))) - Phi_approx((m_min[i] - delay[i,t])/(sqrt(s_min[i]^2+std_dev[i]^2)));
	log_lik[i] += bernoulli_lpmf(choice[i,t] |theta);
	// generate posterior prediction for current trial
	y_pred[i,t] = bernoulli_rng(theta);
      }
    }
  }
}
