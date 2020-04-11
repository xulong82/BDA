# Response rate as an endpoint in cancer trial is a proportion variable 
# How to do statistical analysis of it?
# Does it follow Normal distribution?

# Simulate binomial variable

library(dplyr)

N = 1e3

y = sapply(1:1e3, function(i) {
  rbinom(n = 1e3, size = 1, p = 0.3) %>% sum
})

hist(y)

# E(y) = N * p
mean(y)
N * 0.3

# Var(y) = N * p * (1 - p)
var(y)
N * 0.3 * 0.7

pi = y / 1e3

hist(pi)

# Var(pi) = Var(y) / N^2

var(pi)
var(y) / 1e3 / 1e3

