# Multiple comparison correction issue

library(rstan)
library(dplyr)

setwd("~/Git/stats/bayes/")

# simulate 10 data points each time
# under frequentist statistics, w alpha = 0.05, type i error rate equals alpha

tau = 1; sigma = 1; N = 5

sim1 = lapply(1:1e3, function(i) {
  theta <- rnorm(N, 0, tau)
  y <- rnorm(N, theta, sigma)
  f <- lm(y ~ 1) %>% summary
  p <- f$coefficients[, "Pr(>|t|)"]
  list(y = y, p = p)
})

pval = lapply(sim1, function(x) x$p) %>% unlist
table(pval < 0.05)

# what happens in the bayesian statistics?
which(pval < 0.05)
data <- list(N = N, y = sim1[[19]]$y) 
lm(data$y ~ 1) %>% summary

# flat prior
mod1 <- rstan::stan_model(file = "./lm1.stan", model_name="model 1")

(fit1 <- sampling(mod1, data, warmup = 3e2, iter = 2e3, chains = 3))
plot(fit1)
summary(fit1)$summary

# w/o any prior specification (flat prior), there is no multiplicity control

# normal prior
mod2 <- rstan::stan_model(file = "./lm2.stan", model_name="model 2")
(fit2 <- sampling(mod2, data, warmup = 3e2, iter = 2e3, chains = 3))
plot(fit2)
summary(fit2)$summary
summary(lm(data$y ~ 1))$coefficient

# what is the "right" prior?

x = rnorm(1e6, 0, 1)
sd(x)
y = lapply(1:1e3, function(i) sample(x, size = 1e2, replace = T))
y = sapply(y, mean)
sd(y)

# under the null hypothesis, p-value follow uniform distribution
f = lapply(1:1e3, function(i) {
  y = sample(x, size = 1e1, replace = T)
  summary(lm(y ~ 1))$coefficient
})

p = sapply(f, function(f1) f1[1, "Pr(>|t|)"])

