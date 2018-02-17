library(invgamma)

x = rinvgamma(1e4, shape = 10, scale = 1)
dx = density(x)
dx$x[which.max(dx$y)]

plot(density(x), xlim = c(0, 1))

1 / (x + 1) = sigma

x = 1 / sigma - 1 
