
x1 = sample(c(0, 1), 100, replace = T)
x2 = rnorm(100, 1, 1)
x3 = sample(c(0, 1), 100, replace = T)

y = rnorm(100, 1, 1)

summary(lm(y ~ x1 * x3))
summary(lm(y ~ x1 * x2 + x1 * x3))
