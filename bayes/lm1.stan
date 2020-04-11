data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha, sigma);
}
