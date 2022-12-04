data {
  int<lower=0> J;
  vector[J] mA;
}
parameters {
  vector[J] mu;
  vector<lower=0>[J] sigma;
}
model {
  // priors
  for (j in 1:J){
    mu[j] ~ normal(0, 100);
    sigma[j] ~ inv_chi_square(0.01);
  }
  // likelihood
  for (j in 1:J)
  mA[j] ~ normal(mu[j], sigma[j]);
}
generated quantities {
  vector[J] ypred;
  vector[J] log_lik;
  for(j in 1:J) {
    log_lik[j] = normal_lpdf(mA[j] | mu[j],sigma[j]);
  }
  
  for (i in 1:12) {
    ypred[i] = normal_rng(mu[i], sigma[i]);
  }
}
