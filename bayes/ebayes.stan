data {
  int<lower=1> N;
  vector[N] y;
  vector[N] g;
  real<lower=0> nu;
  real<lower=0> tau;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ scaled_inv_chi_square(nu, tau);
  y ~ normal(alpha + beta * g, sigma);
}

