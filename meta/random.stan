data {
  int<lower=0> J; // number of studies 
  real y[J]; // estimated effect sizes 
  real<lower=0> sigma[J]; // se of estimated effect sizes 
}
parameters {
  real mu; 
  real<lower=0> tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}

