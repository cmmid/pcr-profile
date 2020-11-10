data {
  int<lower=0> N; // number of data points
  int<lower=0> P; // number of patients
  int patient_ID[N];
  vector[N] day_of_test; // day of test (integer)
  int test_result[N]; // boolean test result
  vector[P] time_first_symptom; // day of first symptomatic test
  vector[P] time_last_asym; // day of last test with no symptoms (currently or previously)
  vector[P] te_upper_bound; // maximum time at which infection must have occured before
  real<lower = 0> lmean; // incubation period parameters
  real<lower = 0> lsd;
}

parameters {
  vector[3] beta;
  vector <lower = 0, upper = 1> [P] prop_Te; 
  real <lower = 0> cutpoint;
}

transformed parameters {
  vector <lower = 0> [P] T_e = (te_upper_bound  .* prop_Te);
  vector[N] diff = day_of_test - T_e[patient_ID];
  vector [N] lp;
  
  for(i in 1:N) {
   lp[i] = bernoulli_logit_lpmf(test_result[i] | beta[1] + beta[2] * (diff[i] - cutpoint) + (diff[i] - cutpoint) * beta[2] * beta[3] * step(diff[i] - cutpoint)) +
   bernoulli_logit_lpmf(0 | beta[1] - cutpoint * beta[2]);
  }
}

model {
  
  prop_Te ~ uniform(0, 1);
  
  // Infection time likelihood
  for(j in 1:P) {
    target += log(lognormal_cdf(time_first_symptom[j] - T_e[j], lmean, lsd) - 
      (time_last_asym[j] - T_e[j] <= 0 ? 0 : lognormal_cdf(time_last_asym[j] - T_e[j], lmean, lsd)));
  }

  // PCR positivity likelihood
  // These priors are tight because there is a level of correlation between beta2 and beta3 that
  // makes it a very difficult posterior for stan to sample from.
  // If the equation is 
  // beta[2] * (i - cutpoint) + (i - cutpoint) * beta[3] * step(i - cutpoint))
  // then the relationship between beta2 and beta3 is a ridge as shown here:
  // https://mc-stan.org/docs/2_24/stan-users-guide/collinearity-section.html#collinearity.section
  // In this current formulation the multiplication causes a banana-shaped relationship between beta2
  // and beta3 which are also notoriously difficult to sample from. This is compounded by the fact
  // that there is a logit transform. To prevent divergent transitions the prior is tight around 0
  // which produces the same results but without divergent transitions. You can verify this
  // yourself by making the priors for beta wider.
  beta ~ std_normal();
  cutpoint ~ normal(5, 5) T[0, ];
  
  target += sum(lp);
}

generated quantities {
  vector[301] p;
  real k;
  for(j in 1:301) {
    k = (j * 0.1) - 0.1;
    p[j] = inv_logit(beta[1] + beta[2] * (k - cutpoint) + (k - cutpoint) * beta[3] * beta[2] * step(k - cutpoint));
  }
}
