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
