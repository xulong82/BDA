data {
  int<lower=1> n;
  vector[n] y;
}

parameters {
  real<lower=2> nu;
  real<lower=0.1> tau;
}

model {
  y ~ scaled_inv_chi_square(nu, tau);
}

