
x = rnorm(1e3, 10, 1)
hist(x)

y = rnorm(1e3, 20, 5)
hist(y)

hist(x + y, breaks = 1e2)
hist(x / y, breaks = 1e2)

?rexp

y = rexp(1e3, 2)

mean(y)
var(y)

n = 1e2
lambda = 5

y1 = rnorm(1e3, mean = 1/lambda, sd = sqrt(1/n/lambda^2))

mean(log(y1))
-log(lambda)

var(log(y1))
1/n

