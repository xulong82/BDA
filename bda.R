# Confidence interval illustration
# Confidence interval gives both p-value and effect size
# CLT: sample mean follows normal distribution

X = rnorm(1e3, 0, 1)
hist(X)

N = 30
(Q = qnorm(1 - 0.05/2))

plot(0, 0, xlim = c(-1, 1), ylim = c(1, 100), xlab = "value", ylab = "index")
for(i in 1:100) {
  x = sample(X, N)
  se = sd(x) / sqrt(N)
  interval = c(mean(x) - Q*se, mean(x) + Q*se)
  mycol = ifelse(min(interval) < 0 & max(interval) > 0, 1, 2)
  lines(interval, c(i, i), col = mycol)
}

# Example 1.4

Pp = 0.5 * 0.5^3 
Py = 0.5 * 0.5^3 + 0.5 * 1

Pp = 0.2 * 1
Py = 0.2 * 1 + 0.8 * 0 

Pp = 3.12e-7 * 0.975
Py = 3.12e-7 * 0.975 + 7.6e-5 * 0.00193 + 6.05e-6 * 0.000143
Pp / Py

x = seq(-5, 5)
y_logit <- exp(x) / (1 + exp(x))
y_probit <- pnorm(x)

cor(y_logit, y_probit)

plot(x, y_logit, type = "l")
lines(x, y_probit, col = "red")

plot(y_logit, x, type = "l")
lines(y_probit, x, col = "red")

# Inverse transform sampling 
# in essence, partition areas under the curve, use cdf to find the maximal x

x = seq(0, 5, by = 0.01)
cdf = 1 - exp(-x) # Exponential distribution
plot(x, cdf, type = "l")

v = sapply(1:1e3, function(x) {
  u = runif(1)
  -log(1 - u)
})

hist(v, n = 100)
hist(rnorm(1e3), n = 100, xlim = c(-5, 5))
hist(rt(1e3, 7), n = 100, xlim = c(-5, 5))

# Andrew's World Cup example

# Binomial distribution
# Probablity density function vs probability distribution function

dbinom(x = 1:30, size = 100, prob = 0.3)
sum(dbinom(x = 1:30, size = 100, prob = 0.3))
plot(1:100, dbinom(x = 1:100, size = 100, prob = 0.3), type = "b")

pbinom(q = 30, size = 100, prob = 0.3)
hist(rbinom(1e3, 1e3, 0.3))

library(gtools)
combinations(10, 2)
permutations(10, 2)

dbinom(2, 10, prob = 0.3)
nrow(combinations(10, 2)) * 0.3^2 * 0.7^8

# Beta distribution
seq = seq(0, 1, by = 1e-2)
plot(seq, dbeta(x = seq, shape1 = 5, shape2 = 1), type = "b")
plot(seq, dbeta(x = seq, shape1 = 10, shape2 = 2), type = "b")
plot(seq, dbeta(x = seq, shape1 = 438, shape2 = 544), type = "b")
plot(seq, dnorm(x = seq, mean = 0.446, sd = 0.016), type = "b")

# Jeffery's invariance argument
# f = 1 / p / (1-p)
# q = h(p) = p^2
# f(q) =  1 / p^2 / (1 - p^2)

# dp / dq = 0.5*q + 1
# 1 / p^0.5 / (1-p^0.5)

# Logistic transformation
x = seq(0, 1, 0.01)
y = log(x / (1-x))
plot(x, y)

# sampling with grid
x = c(-0.86, -0.30, -0.05, 0.73)
n = c(5, 5, 5, 5)
y = c(0, 1, 3, 5)

log_post <- function (a,b,y,n,x) {
  loglik = sum (dbinom (y, n, invlogit(a+b*x), log=TRUE))
  return(exp(loglik))
}

a0 = seq(-5, 10, 0.1)
b0 = seq(-10, 40, 0.1)
grid.points <- expand.grid(a0, b0)

post.ord <- apply(grid.points, 1, function(t) log_post(a = t[1], b = t[2], y = y, n = n, x = x))

sample.indices <- sample(1:nrow(grid.points), size = 10000, replace = T, prob = post.ord)
sim.posterior <- grid.points[sample.indices, ]
hist(sim.posterior[, 1])
hist(sim.posterior[, 2])

# contour plot is non-trivial
contour(x = grid.points[, 1], y = grid.points[, 2], z = post.ord)

triangle.prior <- function(x) {
  if (x >= 0 && x < 0.25)
    8 * x
  else if (x >= 0.25 && x <= 1)
    8/3 - 8 * x/3
  else 0
}

posterior.function <- function(theta, n, y) {
  (theta^y) * (1 - theta)^(n - y) * triangle.prior(theta)
}

x = seq(-1, 2, 0.01)
y = sapply(x, triangle.prior)
plot(x, y)

m <- 100
grid.points <- seq(from = 0, to = 1, length.out = m)

unnormal.post.ord <- posterior.function(theta = grid.points, n = 500, y = 285)

k <- 1/m
normal.constant <- sum(k * unnormal.post.ord)
post.ord <- unnormal.post.ord/normal.constant
plot(grid.points, post.ord, type = "l", col = "red")

set.seed(12345)
posterior.triangle.1 <- sample(grid.points, size = 10000, replace = T, prob = post.ord)
hist(posterior.triangle.1, xlim = c(0, 1))

poisson.posterior <- function(theta, y, x, prior.mean.a, prior.var.a, prior.mean.b, prior.var.b) {
  a <- theta[1]
  b <- theta[2]
  lambda <- exp(a + b * x)
  log.like <- sum(dpois(y, lambda = lambda, log = T))
  log.prior.a <- dnorm(a, mean = prior.mean.a, sd = sqrt(prior.var.a), log = T)
  log.prior.b <- dnorm(b, mean = prior.mean.b, sd = sqrt(prior.var.b), log = T)
  log.post <- log.like + log.prior.a + log.prior.b
  return(exp(log.post))
}

y.vec <- sanction$num
x.vec <- sanction$coop
mu.a <- mu.b <- 0
sigma2.a <- sigma2.b <- 20
mle <- glm(num ~ coop, data = sanction, family = poisson)$coef
mle.se <- summary(glm(num ~ coop, data = sanction, family = poisson))$coef[, 2]

grid.a <- seq(from = mle[1] - 5 * mle.se[1], to = mle[1] + 5 * mle.se[1], length.out = 200)
grid.b <- seq(from = mle[2] - 5 * mle.se[2], to = mle[2] + 5 * mle.se[2], length.out = 200)
grid.points <- expand.grid(grid.a, grid.b)

post.ord <- apply(grid.points, MARGIN = 1, FUN = poisson.posterior,
                  y = y.vec, x = x.vec, prior.mean.a = mu.a, prior.var.a = sigma2.a,
                  prior.mean.b = mu.b, prior.var.b = sigma2.b)

sample.indices <- sample(1:nrow(grid.points), size = 10000, replace = T, prob = post.ord)
sim.posterior <- grid.points[sample.indices, ]
