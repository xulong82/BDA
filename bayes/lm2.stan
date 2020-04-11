data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 1);
  y ~ normal(alpha, sigma);
}
