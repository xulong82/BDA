# Exponential distribution model the distribution of the time until next event occurs

# Formal definition of exponential distribution is: 
# the probability distribution of the time *between* the events in a Poisson process.

# T ~ Exp(lambda)
# lambda: event rate, or number of events per unit time

# The exponential distribution, therefore, has the somewhat remarkable property that: 
# we arrive at the exact same inference if we follow d subjects until all have failed or 
# if we follow some larger number n until d have failed


n = 3
lambda = 3

t1 = rexp(1e3, rate = n * lambda)
hist(t1, breaks = 1e2)

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

# delta method for variance-stabilizing

x1 = rexp(n = 1e4, rate = 2)
mean(x1)
var(x1)

xh = sapply(1:1e3, function(x) mean(rexp(n = 10, rate = 2)))
mean(xh)
var(xh)

xh2 = sapply(1:1e4, function(x) log(mean(rexp(n = 10, rate = 2))))
mean(xh2)
-log(2)
var(xh2)
