data {
  vector[2] shape;
}

parameters {
  real<lower=0, upper=1> y;
}

model {
# y ~ beta(shape[1], shape[2]);
  target += beta_lpdf(y | shape[1], shape[2]);
}

