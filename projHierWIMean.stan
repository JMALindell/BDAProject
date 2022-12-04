data {
  int<lower=0> J;
  vector[J] mA;
  int<lower=1,upper=J> gI[J];
}
parameters {
  real mu0; //prior mean
  real<lower=0> sigma0;
  vector[J] mu;
  real<lower=0> sigma;
}
model {
  // priors
  mu0 ~ normal(0,100);
  sigma0 ~ inv_chi_square(0.01);
  sigma ~ inv_chi_square(0.01);  
  for (j in 1:J){
    mu[j] ~ normal(mu0,sigma0);
  }
  // likelihood
  mA ~ normal(mu[gI],sigma);
}
generated quantities {
  vector[J+1] ypred;
  
  vector[J] log_lik;
  for(j in 1:J) {
    log_lik[j] = normal_lpdf(mA[j] | mu[j],sigma);
  }
  
  for (i in 1:J) {
    ypred[i] = normal_rng(mu[i],sigma);
  }
  ypred[13] = normal_rng(mu0,sigma0);
}
