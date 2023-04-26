data {
  int<lower=1> T;
  real std_dev;
  int<lower=0, upper=1> choice[T];
  real delay[T];
}


parameters {
  real m_min;
  real m_plus;
  real<lower = 0> s_min;
  real<lower = 0> s_plus;
}


model {
  m_min  ~ normal(0, 10);
  m_plus  ~ normal(0,10);
  s_min  ~ normal(0, 1);
  s_plus  ~ normal(0,1);

  for (t in 1:T) {
    real theta;
    theta =   Phi_approx((m_plus - delay[t])/(sqrt(s_plus^2+std_dev^2))) - Phi_approx((m_min - delay[t])/(sqrt(s_min^2+std_dev^2)));
    choice[t] ~ bernoulli(theta);
  }
}


generated quantities {
  // For log likelihood calculation
  real log_lik;

  // For posterior predictive check
  real y_pred[T];

  // Set all posterior predictions to -1 (avoids NULL values)
  for (t in 1:T) {
    y_pred[t] = -1;
  }
  
  { // local section, this saves time and space
    real theta;

    // Initialize values
    log_lik = 0;
    for (t in 1:T) {
      theta =   Phi_approx((m_plus - delay[t])/(sqrt(s_plus^2+std_dev^2))) - Phi_approx((m_min - delay[t])/(sqrt(s_min^2+std_dev^2)));
      log_lik += bernoulli_lpmf(choice[t] |theta);
       // generate posterior prediction for current trial
      y_pred[t] = bernoulli_rng(theta);
     }
  }
}

