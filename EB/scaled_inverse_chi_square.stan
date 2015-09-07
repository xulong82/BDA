data {
  real<lower=0> nu;
  real<lower=0> tau;
}

parameters {
  real<lower=0> sigma;
}

model {
  sigma ~ scaled_inv_chi_square(nu, tau);
}

