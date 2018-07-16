data {
  int<lower=0> J; // number of studies 
  real y[J]; // estimated effect sizes 
  real<lower=0> sigma[J]; // se of estimated effect sizes 
}
parameters {
  real mu; 
}
model {
  y ~ normal(mu, sigma);
}

