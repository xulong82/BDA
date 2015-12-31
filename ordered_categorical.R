stan_ordered <- "

data {
  int<lower=1> N;  // Sample number
  int<lower=1> K;  // Ad categories
  int<lower=1> D;  // Predictors
  row_vector[D] cov[N]; // Covariates
  int<lower=1,upper=K> Ad[N];  // Ad status
}

parameters {
  ordered[K-1] c;
  vector[D] beta;
} 

model {
  c ~ cauchy(3, 1);
  beta ~ cauchy(0, 1);

  for (n in 1:N)
    Ad[n] ~ ordered_logistic(cov[n] * beta, c);
}

"

stan_ordered <- stan_model(model_code = stan_ordered)

(fit <- sampling(stan_ordered, data = dat, chain = 3, iter = 400, warmup = 200))
plot(fit)

(fit <- sampling(stan_ordered, data = dat, chain = 3, 
                 init = list(c1 = init, c2 = init, c3 = init),
                 iter = 400, warmup = 200, show_messages = F)) # 15 sec 

plot(fit)

samples = extract(fit) # samples

x = summary(fit)$summary[, "mean"]
init = as.list(x)
init$lp__ = NULL

load("~/Dropbox/GitHub/Adsp/data/mdata.rdt")
cov <- mdata[c("Age", "Sex", "Apoe2", "Apoe4")]
cov$Sex = as.integer(cov$Sex) - 1
dat = list(N = 576, K = 4, D = 4, cov = cov, Ad = as.numeric(mdata$AD2))
