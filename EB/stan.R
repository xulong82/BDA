library(rstan)

setwd("~/Dropbox/Stan")

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
		    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = '8schools.stan', data = schools_dat, iter = 1000, chains = 4)

insurance = read.delim("slr06.txt")

model = stan_model(file = "insurance.stan")
dat = list(N = nrow(insurance), X = insurance$X, Y = insurance$Y)
fit = sampling(model, data = dat, chains = 2)

print(fit)

# residual variance: sqrt(sigma)

lm(Y ~ X, data = insurance) %>% summary
lm(Y ~ 1, data = insurance) %>% summary
