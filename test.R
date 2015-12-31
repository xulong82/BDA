stan_glm <- "

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  int<lower=1> D;  // Predictors
  row_vector[D] cov[N];  // Covariates
  cov_matrix[N] Sigma;  // Kinship
  int<lower=1,upper=K> Ad[N];  // Ad status
}

transformed data {
  matrix[N, N] L;
  L <- cholesky_decompose(Sigma);
}

parameters {
  ordered[K-1] c;
  vector[D] beta;
  vector[N] z; 
} 

transformed parameters {
  vector[N] random;
  random <- L * z;
}

model {
  c ~ cauchy(3, 1);
  beta ~ cauchy(0, 1);
  z ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(cov[n] * beta + random[n], c);
}

"

stan_glm <- stan_model(model_code = stan_glm)

(fit <- sampling(stan_glm, data = dat, chain = 3, iter = 400, warmup = 200))
plot(fit)

(fit <- sampling(stan_ordered, data = dat, chain = 3, 
                 init = list(c1 = init, c2 = init, c3 = init),
                 iter = 400, warmup = 200, show_messages = F)) # 15 sec 
plot(fit)

