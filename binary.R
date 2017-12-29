x = c(rep("A", 10), rep("B", "15"), rep("C", 20))

x = sample(x)

x = sapply(c("A", "B", "C"), function(xx) as.numeric(x == xx))

summary(lm(x[, 1] ~ x[, 2] + x[, 3]))

y = rnorm(45, mean = 1, sd = 2)
