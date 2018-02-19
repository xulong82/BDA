library(lmtest)

# there is no error terms in logistic regression models, because there is no error terms in binomial distribution model
# but there are still residuals, variances, and such on

# Frank Harrell
# There is no error term in the Bernoulli distribution, there's just an unknown probability. The logistic model is a probability model.

# deviance, residual variance

x = rnorm(100, 0, 1)
y = 1 + 2 * x + rnorm(100, 0, 1)
plot(x, y)

z = glm(y ~x)

sd(z$residuals)
sd(z$residuals)^2
var(z$residuals)

z$deviance
sum(z$residuals^2)

yp = predict(z)
res = abs(y - yp)
sum(res^2)

z$deviance / (100 - 1)
var(z$residuals)

# deviance and LRT relationship

z1 = glm(y ~ 1)
z2 = glm(y ~ 1 + x)

lrtest(z1, z2)
z1$deviance - z2$deviance

# does not hold in gaussian model

# logistic model

x = rbinom(n = 100, size = 1, prob = 0.8)

z1 = glm(x ~ 1, family = "binomial")
z2 = glm(x ~ rnorm(100, 1, 1), family = "binomial")

z1$deviance - z2$deviance
lrtest(z2, z1)

z2$deviance
-2 * logLik(z2)

# hold in logistic model

z1 = glm(x ~ 1, family = "gaussian")
z2 = glm(x ~ rnorm(100, 1, 1), family = "gaussian")

z1$deviance - z2$deviance
lrtest(z2, z1)

z2$deviance
z2$deviance / 0.1438492

-2 * logLik(z2)

# again, does not hold in gaussian model

var(z2$residuals)
var(x)

z2$fitted.values
z2$linear.predictors

# turns out there are scaled and unscaled deviance
# it is the scaled deviance that equates the lrt, but not necessarily the unscaled deviance
# scaled deviance is calcualted by dividing unscaled deviance by a constant parameter: the dispersion parameter
# this dispersion parameter is 1 for binomial and possion distribution, so the scaled and unscaled deviance are the same
# but, the dispersion parameter is not 1 for gaussian, and this is the reason
