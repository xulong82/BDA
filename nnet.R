library(nnet)

n <- 1000
df1 <- data.frame(x1=runif(n,0,100), x2=runif(n,0,100))
df1 <- transform(df1, y=1+ifelse(100 - x1 - x2 + rnorm(n, sd=10) < 0, 0, ifelse(100 - 2*x2 + rnorm(n,sd=10) < 0, 1, 2)), set="Original")

model <- multinom(y ~ x1 + x2, df1)

predict(model, df1, "probs")

x = rnorm(10, 0, 1)
sort(x)

ecdf(x)
plot(ecdf(x))
