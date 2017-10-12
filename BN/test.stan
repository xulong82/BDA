data {
  int<lower=2> N; # sample number 
  int<lower=2> K; # variable number
  matrix[K, N] Y; # data matrix
  matrix[K, K] X; # design matrix 
}

parameters {
  vector[K] alpha;
  vector[K] sigma;
  matrix[K, K] beta0;
}

transformed parameters {
  matrix[K, K] beta = X .* beta0;
}

model {
  real ps[K];

  alpha ~ normal(0, 10);
  sigma ~ inv_gamma(.001, .001);

  for (n in 1:K) 
    beta0[, n] ~ normal(0, 1); 

  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = normal_lpdf(Y[k, n] | alpha[k] + beta[1, k] * Y[1, n] + beta[2, k] * Y[2, n] + beta[3, k] * Y[3, n] + beta[4, k] * Y[4, n], sigma[k]);
    }
    target += log_sum_exp(ps);
  }

}

generated quantities {
  matrix[K, N] lp;

  for (n in 1:N) {
    for (k in 1:K) {
      lp[k, n] = normal_lpdf(Y[k, n] | alpha[k] + beta[1, k] * Y[1, n] + beta[2, k] * Y[2, n] + beta[3, k] * Y[3, n] + beta[4, k] * Y[4, n], sigma[k]);
    }
  }
}

