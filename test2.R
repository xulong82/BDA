# Let x ~ Exp(lambda)
# x_1 ~ Exp(n * lambda)

x1 = rexp(n = 1e3, rate = 1)
x2 = rexp(n = 1e3, rate = 10)

hist(x1)
hist(x2)

mean(x1)
mean(x2)

t1 = sapply(1:1e3, function(x) min(rexp(n = 10, rate = 1)))
hist(t1)
mean(t1)
