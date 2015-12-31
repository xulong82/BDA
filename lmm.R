stan_glm <- "

data {
  int<lower=1> N;  // Sample number
  int<lower=1> D;  // Predictors
  row_vector[D] cov[N];  // Covariates
  cov_matrix[N] Sigma;  // Kinship
  vector[N] Ad;  // Ad status
}

transformed data {
  matrix[N, N] L;
  L <- cholesky_decompose(Sigma);
}

parameters {
  real alpha;
  vector[D] beta;
  real<lower=machine_precision()> theta_e;
  real<lower=machine_precision()> theta_u;
  vector[N] u; 
} 

model {
  alpha ~ cauchy(0, 1);
  beta ~ cauchy(0, 1);
  u ~ multi_normal_cholesky(rep_vector(0, N), theta_u * L);

  {
    vector[N] mu;
    
    for (n in 1:N)
      mu[n] <- alpha + cov[n] * beta + u[n];

    Ad ~ normal(mu, theta_e);
  }

}

"

stan_glm <- stan_model(model_code = stan_glm)

(fit <- sampling(stan_glm, data = dat, chain = 3, iter = 400, warmup = 200))
plot(fit)

(fit <- sampling(stan_glm, data = dat, chain = 3, 
                 init = list(c1 = init, c2 = init, c3 = init),
                 iter = 400, warmup = 200, show_messages = F))
plot(fit)

samples = extract(fit) # samples

z = samples$z
random = samples$random

var(random)[1:10, 1:10]
var(z)[1:10, 1:10]

x = summary(fit)$summary[, "mean"]
init = as.list(x)
init$lp__ = NULL

load("~/Dropbox/GitHub/Adsp/data/mdata.rdt")
load("~/Dropbox/GitHub/Adsp/data/kinship.rdt")
cov <- mdata[c("Age", "Sex", "Apoe2", "Apoe4")]
cov$Sex = as.integer(cov$Sex) - 1
Sigma <- kinship$autosome
Sigma[Sigma < 0] <- 0

dat = list(N = 576, D = 4, cov = cov, Sigma = Sigma, Ad = mdata$AD1)
