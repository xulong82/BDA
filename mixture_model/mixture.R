# mixture model: model multiple data generating processes
# http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html

library(rstan)
library(parallel)

rm(list = ls())
setwd("~/Projects/stats")

# a mixture of 2 gaussian

mu <- c(-2.75, 2.75)
sigma <- c(1, 1)
lambda <- 0.4

# simulation

set.seed(689934)

N <- 1e3
z <- rbinom(N, 1, lambda) + 1
y <- rnorm(N, mu[z], sigma[z])

# simulate beta samples

idata = list(shape = c(5, 5))
fit <- stan(file = "beta.stan", data = idata)
