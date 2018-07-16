data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Category number
  int<lower=1> D;  // Predictor number
  row_vector[D] cov[N];  // Covariates
  matrix[N, N] Sigma;  // Kinship
  int<lower=1,upper=K> Ad[N];  // Response
}

parameters {
  ordered[K-1] c;
  vector[D] beta;
  vector[N] z; 
} 

transformed parameters {
  vector[N] random;
  random <- Sigma * z;
}

model {
  c ~ normal(3, 1);
  beta ~ normal(0, 1);
  z ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(cov[n] * beta + random[n], c);
}

# real<lower=machine_precision()> scale;
# random <- scale * Sigma * z;

